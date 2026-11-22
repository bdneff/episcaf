"""Check uploaded coordinates against the fixed pilot reference."""
from episcaf_service import lysozyme_job


def matches_reference(structure, route, chain):
    if chain != "Y":
        return False
    reference = lysozyme_job.REFERENCE.read_text() if route == "antibody" else lysozyme_job.reference_pdb()
    chains = "HLY" if route == "antibody" else "Y"
    def atoms(text):
        return [line.rstrip() for line in text.splitlines() if line.startswith("ATOM  ") and line[21:22] in chains]
    return atoms(structure) == atoms(reference)
