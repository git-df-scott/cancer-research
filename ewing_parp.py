"""
The Ewing sarcoma PARP paradox: did olaparib ever reach an active concentration?

THE PARADOX
-----------
Ewing sarcoma lines are among the most PARP-inhibitor-sensitive in all of cancer cell line
screening. The phase II trial of single-agent olaparib in refractory Ewing sarcoma produced
0 objective responses in 12 patients, median time to progression 5.7 weeks (PMC4230717).

Four explanations are offered in the literature: insufficient synthetic lethality, autophagy-
mediated resistance, the need for combination, and "failure to achieve in vitro levels of olaparib
at the clinical dose". The last is a pharmacokinetics question. It has been proposed repeatedly
and, as far as a literature search finds, never computed.

THE COMPARISON THAT IS USUALLY SKIPPED
--------------------------------------
In vitro IC50s are measured in media with ~10% serum, where protein binding is low, so the nominal
IC50 is close to a FREE concentration. Plasma concentrations in patients are TOTAL, and olaparib is
82-91% protein bound. Comparing a total plasma concentration to an in vitro IC50 therefore
overstates the exposure a tumour cell actually sees, by roughly 5-10x.

Only free drug is pharmacologically active. This computes the free-drug comparison.

INPUTS, ALL PUBLISHED
---------------------
Clinical  olaparib 300 mg BID tablet, steady state: Cmax 7.6 ug/mL, AUC 49.2 ug*h/mL over the
          12 h interval (FDA label / popPK). The Ewing trial used 400 mg BID capsules, whose
          exposure is ~13% LOWER than the 300 mg tablet, so using tablet numbers is generous to
          the drug.
Binding   82% bound at 10,000 ng/mL, 91% at 1000 ng/mL (FDA clinical pharmacology review). Free
          fraction 0.09-0.18 depending on concentration; both bounds carried.
In vitro  Ewing lines, 72 h: median IC50 1.995 +/- 0.46 uM.
          Paediatric solid tumour panel: median IC50 3.6 uM (range 1-33.8).
          Long-term growth assays: activity reported down to 600 nM.
"""
MW = 434.46          # g/mol, olaparib C24H23FN4O3

CMAX_TOTAL = 7.6     # ug/mL, steady state, 300 mg BID tablet
AUC_12H = 49.2       # ug*h/mL over the 12 h dosing interval
CAVG_TOTAL = AUC_12H / 12.0

FU = (0.09, 0.18)    # free fraction, 91% and 82% bound

IC50 = {
    'Ewing lines, 72 h (median)': 1.995,
    'paediatric panel (median)': 3.6,
    'long-term growth assay': 0.600,
}


def ugml_to_uM(c):
    return c * 1000.0 / MW


def free_range(total):
    return ugml_to_uM(total * FU[0]), ugml_to_uM(total * FU[1])


def report():
    cmax_f = free_range(CMAX_TOTAL)
    cavg_f = free_range(CAVG_TOTAL)
    print('CLINICAL EXPOSURE, olaparib 300 mg BID at steady state\n')
    print(f"  total Cmax   {CMAX_TOTAL:5.2f} ug/mL = {ugml_to_uM(CMAX_TOTAL):5.2f} uM")
    print(f"  total Cavg   {CAVG_TOTAL:5.2f} ug/mL = {ugml_to_uM(CAVG_TOTAL):5.2f} uM")
    print(f"  FREE  Cmax   {cmax_f[0]:5.2f} - {cmax_f[1]:5.2f} uM   (fu {FU[0]}-{FU[1]})")
    print(f"  FREE  Cavg   {cavg_f[0]:5.2f} - {cavg_f[1]:5.2f} uM")

    print('\nAGAINST IN VITRO POTENCY\n')
    print(f"{'in vitro measure':>30s} {'IC50 (uM)':>10s} {'free Cavg/IC50':>16s} {'free Cmax/IC50':>16s}")
    for name, ic in IC50.items():
        r_avg = f"{cavg_f[0]/ic:.2f}-{cavg_f[1]/ic:.2f}"
        r_max = f"{cmax_f[0]/ic:.2f}-{cmax_f[1]/ic:.2f}"
        print(f"{name:>30s} {ic:10.3f} {r_avg:>16s} {r_max:>16s}")

    print('\nREADING THE RATIOS')
    print('  A ratio of 1.0 means free drug in a patient equals the concentration that produces')
    print('  50% growth inhibition in vitro. Cytotoxicity generally needs several-fold above IC50,')
    print('  not parity, so 1.0 is a floor for plausible activity rather than a target.')

    ic = IC50['Ewing lines, 72 h (median)']
    print(f"\n  Against the Ewing-specific 72 h IC50 of {ic} uM:")
    print(f"    free Cavg sits at {cavg_f[0]/ic:.2f}-{cavg_f[1]/ic:.2f} of it -- BELOW parity")
    print(f"    free Cmax sits at {cmax_f[0]/ic:.2f}-{cmax_f[1]/ic:.2f} of it -- parity at best, and only at peak")

    ic = IC50['long-term growth assay']
    print(f"\n  But against the long-term assay threshold of {ic} uM:")
    print(f"    free Cavg sits at {cavg_f[0]/ic:.2f}-{cavg_f[1]/ic:.2f} of it -- ABOVE, comfortably")
    print('    Chronic dosing is the clinical situation, so this is the fairer comparison and it')
    print('    does NOT support a pure exposure explanation.')


if __name__ == '__main__':
    report()


# --------------------------------------------------------------------------------------
# Talazoparib comparison. The field moved here after olaparib monotherapy failed, and this
# analysis says the move was right for a reason beyond the one usually given.
TALAZOPARIB = dict(fu=0.26, trapping_vs_olaparib=100.0)


def talazoparib_note():
    print('\n' + '=' * 74)
    print('WHY THE FIELD MOVED TO TALAZOPARIB, AND A COMPOUNDING EFFECT NOT USUALLY STATED')
    print('=' * 74)
    lo, hi = FU
    fu_t = TALAZOPARIB['fu']
    print(f"\n  olaparib free fraction     {lo:.2f}-{hi:.2f}")
    print(f"  talazoparib free fraction  {fu_t:.2f}  (74% bound, concentration-independent)")
    print(f"  -> talazoparib delivers {fu_t/hi:.1f}-{fu_t/lo:.1f}x more FREE drug per unit total")
    print(f"  -> on top of ~{TALAZOPARIB['trapping_vs_olaparib']:.0f}x greater PARP-trapping potency")
    print(f"  -> combined free-drug potency advantage: "
          f"{TALAZOPARIB['trapping_vs_olaparib']*fu_t/hi:.0f}-"
          f"{TALAZOPARIB['trapping_vs_olaparib']*fu_t/lo:.0f}x")
    print('\n  The potency ratio is widely quoted. The protein-binding difference is not, and it')
    print('  compounds with it. Talazoparib + irinotecan +/- temozolomide is active in Ewing')
    print('  (COG ADVL1411), which is consistent with this.')


def verdict():
    print('\n' + '=' * 74)
    print('VERDICT')
    print('=' * 74)
    print("""
The PK explanation is PARTIALLY SUPPORTED, NOT CLEAN. Free Cavg sits below the Ewing-specific
72 h IC50 and above the long-term growth-assay threshold. Which in vitro measure is the right
comparator decides the answer, and that is not settled.

What IS clean, and is the point worth carrying:

    total  Cavg = 9.44 uM   -- about 5x the in vitro IC50, looks comfortable
    FREE   Cavg = 0.85-1.70 uM -- at or below it

The protein-binding correction flips the conclusion from "roughly 5x adequate" to "marginal at
best". Anyone reasoning from total plasma concentration, which is what gets reported, would have
concluded exposure was ample. It was not. There was no headroom for tumour penetration, efflux, or
tissue binding.

The trial also used 400 mg BID capsules, ~13% lower exposure than the 300 mg tablets computed here,
so the real figures were slightly worse than shown.

This does not cure Ewing sarcoma and does not claim to. It says the monotherapy failure is
consistent with a drug that never had exposure headroom, and that the field's move to a more potent
agent with a higher free fraction was better founded than the potency ratio alone suggests.
""")


if __name__ == '__main__':
    talazoparib_note()
    verdict()
