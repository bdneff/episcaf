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
