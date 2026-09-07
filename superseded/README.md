# Superseded: the first attempt, and why it was wrong

This directory is kept deliberately. It is the failed first attempt at the question, preserved so the
reasoning trail is auditable.

`exp8_follicular.py` bolted a multi-follicle geometry onto the parent project's **solid-tumour** model
(oxygen diffusing from a box boundary, a drug that kills cells attempting division, resistance as an
evolving scalar trait) and asked whether the adaptive-therapy benefit survived in lymphoma-like
architecture. It ran 864 simulations. **Four of its five pre-registered predictions failed.**

Two separate problems, both recorded here rather than quietly dropped:

1. **A design flaw.** The implantation routine placed the resistant clone at the deepest point of
   whichever nodule it chose, so implant depth barely varied across conditions (7.6 to 9.9 sites).
   The prediction that benefit scales with depth was therefore untestable as built. The compact
   tumour also had 20.8% interior vacancies because the box was enlarged from L=120 to L=160,
   starving the core more than in every prior experiment in the parent project.
2. **A wrong abstraction, which matters more.** Follicular lymphoma relapse arises from a
   pre-existing common-progenitor reservoir that diverges from the treated bulk, not from a resistant
   clone evolving inside the lesion. Germinal-centre hypoxia is physiological signalling in a
   ~100 um structure, not diffusion-limited necrosis. The drugs that work are immune cells, not
   diffusing molecules. So even a well-executed version of this experiment would have been answering
   the wrong question.

See `../UNDERSTANDING_NHL.md` sections 4, 5 and 8 for the biology that made this clear, and the
project README for what replaced it.
