// Private endpoints always belong to the site that served this page.
async function pilotRequest(path, options = {}) {
  const response = await fetch("/pilot/" + path, {
    ...options, mode: "same-origin", credentials: "same-origin", redirect: "error", cache: "no-store"
  });
  if (!response.ok) {
    const body = await response.json().catch(() => ({}));
    throw new Error(typeof body.detail === "string" ? body.detail : "The pilot request failed.");
  }
  return response;
}

// Both pages sign in on the current origin. Results also refreshes the run list.
function pilotSession(onRuns, onSignedOut) {
  const form = document.getElementById("signin");
  const signout = document.getElementById("signout");
  const message = document.getElementById("authmsg");
  let revision = 0;
  function signedOut(text) {
    form.hidden = false;
    signout.hidden = true;
    message.textContent = text;
    onSignedOut();
  }
  async function refresh() {
    const current = ++revision;
    try {
      const runs = await (await pilotRequest("runs")).json();
      if (current !== revision) return;
      form.hidden = true;
      signout.hidden = false;
      message.textContent = "";
      await onRuns(runs);
    } catch (error) {
      if (current === revision) signedOut(error.message);
    }
  }
  form.onsubmit = async event => {
    event.preventDefault();
    const input = document.getElementById("token");
    const token = input.value.trim();
    input.value = "";
    try {
      await pilotRequest("session", {method:"POST", headers:{Authorization:"Bearer " + token}});
      await refresh();
    } catch (error) { signedOut(error.message); }
  };
  signout.onclick = async () => {
    try {
      await pilotRequest("session", {method:"DELETE"});
      ++revision;
      signedOut("Signed out.");
    } catch (error) { message.textContent = error.message; }
  };
  return refresh;
}
