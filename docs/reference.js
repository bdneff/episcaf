// Linkable panels with ordinary document navigation when JavaScript is absent.
const sectionNav = document.querySelector("[data-section-nav]");
const sectionLinks = Array.from(sectionNav.querySelectorAll("a"));
const sections = Array.from(document.querySelectorAll("[data-section]"));

function hashTarget(hash) {
  try {
    return document.getElementById(decodeURIComponent(hash.slice(1)));
  } catch {
    return null;
  }
}

function showSection(hash, { focus = false, scroll = true } = {}) {
  const target = hashTarget(hash);
  const requested = target?.closest("[data-section]");
  const section = requested && !requested.hasAttribute("data-section-unavailable")
    ? requested : sections.find(panel => !panel.hasAttribute("data-section-unavailable"));

  for (const panel of sections) {
    panel.hidden = panel !== section;
  }
  for (const link of sectionLinks) {
    if (link.hash.slice(1) === section.getAttribute("aria-labelledby")) {
      link.setAttribute("aria-current", "location");
    } else {
      link.removeAttribute("aria-current");
    }
  }

  for (let disclosure = target?.closest("details"); disclosure;
    disclosure = disclosure.parentElement.closest("details")) {
    disclosure.open = true;
  }

  if (!scroll || !target || (requested && requested !== section)) return;

  if (focus) {
    target.focus({ preventScroll: true });
    target.scrollIntoView({ block: "nearest" });
  } else {
    target.scrollIntoView({ block: "nearest" });
  }
}

function navigateToSection(hash) {
  if (location.hash !== hash) history.pushState(null, "", hash);
  showSection(hash, { focus: true });
}

sectionNav.addEventListener("click", (event) => {
  const link = event.target.closest("a");
  if (!link || event.defaultPrevented || event.button !== 0 ||
    event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) {
    return;
  }
  event.preventDefault();
  navigateToSection(link.hash);
});

window.addEventListener("hashchange", () => showSection(location.hash));
showSection(location.hash);
