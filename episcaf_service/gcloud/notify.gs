// Set REQUEST_BUCKET and OPERATOR_EMAIL in the project's Script Properties.
function notifyOperator() {
  const lock = LockService.getScriptLock();
  if (!lock.tryLock(1000)) return;
  try {
    const config = PropertiesService.getScriptProperties();
    const bucket = config.getProperty('REQUEST_BUCKET');
    const recipient = config.getProperty('OPERATOR_EMAIL');
    if (!bucket || !recipient) throw Error('Set REQUEST_BUCKET and OPERATOR_EMAIL.');
    const base = 'https://storage.googleapis.com/storage/v1/b/' + encodeURIComponent(bucket)
      + '/o?prefix=inbox%2F&fields=nextPageToken,items(name)';
    const options = {headers: {Authorization: 'Bearer ' + ScriptApp.getOAuthToken()}};
    let page = '';
    do {
      const response = UrlFetchApp.fetch(base + '&pageToken=' + encodeURIComponent(page), options);
      const list = JSON.parse(response.getContentText());
      for (const item of list.items || []) {
        const match = /^inbox\/([\w-]+)\/([a-f0-9]{32})\/(reference|review)\.zip$/.exec(item.name);
        if (!match || config.getProperty('sent.' + match[2])) continue;
        const [, owner, id, kind] = match;
        const ref = kind === 'reference';
        const op = 'python3 episcaf_service/gcloud/pilot_transfer.py ';
        const body = [
          owner + ' submitted request ' + id + '.',
          ref ? '3HFM reference: 16 backbone attempts, up to 8 retained.' : 'Structure submitted for lab review.',
          'Open the submitted ZIP with your lab Google account:',
          'https://storage.cloud.google.com/' + encodeURIComponent(bucket) + '/' + item.name,
          'Open request.json inside the ZIP for the filename, chain and selected regions.'
        ];
        if (ref) body.push('To transfer the reference request, run this from your laptop checkout:',
          op + 'fetch ' + owner + ' ' + id,
          'After review, on Gemini:',
          'cd ~/.episcaf-pilot/inbox/' + id + ' && sbatch --parsable ~/episcaf-advisor/episcaf_service/handoff.sbatch',
          'From your laptop, record the numeric job ID:',
          op + 'receipt ' + owner + ' ' + id + ' JOB_ID',
          'When the job finishes: ' + op + 'publish ' + owner + ' ' + id + ' JOB_ID');
        MailApp.sendEmail(recipient, 'episcaf request from ' + owner + ' · ' + id, body.join('\n\n'));
        config.setProperty('sent.' + id, '1');
      }
      page = list.nextPageToken || '';
    } while (page);
  } finally {
    lock.releaseLock();
  }
}
