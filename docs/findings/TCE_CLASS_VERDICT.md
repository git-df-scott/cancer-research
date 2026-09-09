# TCE class-wide PK claim: NOT PUBLISHABLE. Killed by three independent checks.

Elapsed time from hypothesis to retraction: about two hours. Recorded because the failure is more
useful than the claim would have been.

## The claim I was going to make

> No approved half-life-extended T-cell engager can produce a treatment-free interval at its
> labelled schedule. The exhaustion-recovery literature rests on blinatumomab, the one agent in the
> class where stopping the drug removes it.

Mono-exponential screen said 8/8 agents fail a 0.111 trough/peak threshold. It looked clean.

## Why it is dead

**1. The screen was biased 2–3× in my favour.** Mono-exponential decline assumes one terminal
phase. Real TCE PK is multi-compartment, so the within-interval fall includes faster distribution
phases and true trough is lower than the screen gives. Against the validated two-compartment
tarlatamab model: 0.099–0.128 actual versus 0.188–0.420 screened. Tarlatamab straddles the
threshold rather than clearly failing it.

**2. A half-life I used was simply wrong.** I had mosunetuzumab at 6–11 days from apparent early
values. Bender et al. (Clin Transl Sci 2024, PMC11134317) report steady-state terminal half-life
**16.1 days**, with time-dependent clearance transitioning over 16.3 days.

**3. The premise fails empirically, and this is the fatal one.** The threshold argument assumes
agents operate near target saturation at peak, so a concentration fall barely moves occupancy.
Bender **measured** CD20 receptor occupancy for mosunetuzumab:

| condition | Day-21 receptor occupancy |
|---|---|
| no residual obinutuzumab | 55% |
| median residual rituximab (10 µg/mL) | ~5% |
| high baseline obinutuzumab | 0.86% |

Mosunetuzumab's CD20 KD is 10.2 µg/mL — weak — and it competes with residual anti-CD20 therapy.
It is nowhere near saturation. Below saturation, occupancy is approximately linear in
concentration, so the 60% concentration fall between doses **is** a 60% occupancy fall. The holiday
exists for that agent in exactly the terms I claimed it could not.

My argument only bites for agents operating near saturation, and peak occupancy is published for
almost none of them.

## And it is less novel than I thought

Bender et al. already provide PK-simulation-based dose-delay guidance (≤14 d proceed; 15–41 d
repeat last dose; ≥42 d repeat step-up), and cite in vitro data on repeat dosing 7 days apart.
The framing is CRS mitigation rather than exhaustion — but interval analysis for half-life-extended
TCEs is not unoccupied territory.

## What actually survives

Two things, both smaller than a finding:

- **Blinatumomab is a structural outlier.** 2-hour half-life on continuous infusion versus 3.8–22
  days for everything else in the class. Uncontroversial once stated.
- **The Philipp washout critique.** TFI benefit was demonstrated with AMG 562 (~210 h half-life) by
  removing drug at reculture. No clinical manoeuvre reproduces that for a 9-day-half-life molecule.
  This one is still clean and still unaddressed anywhere I can find.

- The EC50-free threshold derivation is correct mathematics with a premise that is empirically
  false for at least one agent. Kept as a tool, not as a result.

## Verdict

**Not publishable, not worth sending to anyone.** A gap statement, not a finding.

The honest summary of the exercise: a hypothesis that survived a prior-art check died within two
hours of contact with the actual pharmacology, killed by reading one paper properly. That is the
process working, and it is also not a result.
