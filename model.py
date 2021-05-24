

def AtoSC_dot(s, t, params):
    ## PARAMETERS
    k_b1, k_d1, k_b2, k_d2, k_b3, k_d3, k_ap, k_ad, k_pt, k_ph, k_dim, k_mon, k_dbnd, k_dunbnd, k_b4, k_d4, k_cat, k_exp, k_lexp, k_deg = params

    ## SPECIES
    C, S, Cp, Sp, pato, Ph, Cpd, CS, CSp, CpS, patoCpd, PhCp, GFP = s

    ## REACTIONS
    # AtoC binds AtoSP: C + Sp -> CSp
    r1 = C * Sp * k_b1
    # AtoC unbinds AtoSP: CSp -> C + Sp
    r2 = CSp * k_d1
    # AtoCP binds AtoS: Cp + S -> CpS
    r3 = S * Cp * k_b2
    # AtoCP unbinds AtoS: CpS -> Cp + S
    r4 = CpS * k_d2
    # AtoC binds AtoS: C + S -> CS
    r5 = C * S * k_b3
    # AtoC unbinds AtoS: CS -> C + S
    r6 = CS * k_d3
    # Acetoacetate phosphorylates AtoS: S -> Sp
    r7 = S * k_ap
    # Dephosphorylation of AtoSP: Sp -> S
    r8 = Sp * k_ad
    # Phosphorylation of AtoC: CSp -> CpS
    r9 = CSp * k_pt
    # Dephosphorylation of AtoC: CpS -> CS (or should this be -> CSp?)
    r10 = CpS * k_ph
    # Dimerisation of AtoCP: Cp + Cp -> Cpd
    r11 = Cp * Cp * k_dim
    # Monomeristaion of AtoCP: Cpd -> Cp + Cp
    r12 = Cpd * k_mon
    # Promoter binding: pato + Cpd -> patoCpd
    r13 = pato * Cpd * k_dbnd
    # Promoter unbinding: patoCpd -> pato + Cpd
    r14 = patoCpd * k_dunbnd
    # Binding of alternative dephosphatase: Cp + Ph -> PhCp
    r15 = Cp * Ph * k_b4
    # Unbinding of alternative dephosphatase: PhCp -> Cp + Ph
    r16 = PhCp * k_d4
    # Alternative dephosphorylation: PhCp -> C + Ph
    r17 = PhCp * k_cat
    # GFP expression: patoCpd -> patoCpd + GFP
    r18 = patoCpd * k_exp
    # Leaky GFP expression: pato -> pato + GFP
    r19 = pato * k_lexp
    # GFP degradation: GFP -> 0
    r20 = GFP * k_deg

    ## ODEs
    dC = r2 + r6 + r17 - r1 - r5
    dS = r4 + r6 + r8 - r3 - r5 - r7
    dCp = r4 + 2 * r12 + r16 - r3 - 2 * r11 - r15
    dSp = r2 + r7 - r1 - r8
    dpato = r14 - r13
    dPh = r16 + r17 - r15
    dCpd = r11 + r14 - r12 - r13
    dCS = r5 + r10 - r6
    dCSp = r1 - r2 - r9
    dCpS = r3 + r9 - r4 - r10
    dpatoCpd = r13 - r14
    dPhCp = r15 - r16 - r17
    dGFP = r18 + r19 - r20

    ds = (dC, dS, dCp, dSp, dpato, dPh, dCpd, dCS, dCSp, dCpS, dpatoCpd, dPhCp, dGFP)

    return ds


def AtoSC_ori_dot(s, t, params):
    ## PARAMETERS
    k_b1, k_d1, k_b2, k_d2, k_b3, k_d3, k_ap, k_ad, k_pt, k_ph, k_b4, k_d4, k_cat, k_deg, k_bnd, k_unbnd, k_nsbnd, k_nsunbnd, k_pmgexp, k_npmgexp, k_mgbnd, k_gexp, k_mat, k_mgdeg = params

    ## SPECIES
    C, S, Cp, Sp, pato, Ph, CS, CSp, CpS, PhCp, patoCp, patoC, GFP, mGFP, R, mGFPR, uGFP = s

    ## REACTIONS
    # AtoC binds AtoSP: C + Sp -> CSp
    r1 = C * Sp * k_b1
    # AtoC unbinds AtoSP: CSp -> C + Sp
    r2 = CSp * k_d1
    # AtoCP binds AtoS: Cp + S -> CpS
    r3 = S * Cp * k_b2
    # AtoCP unbinds AtoS: CpS -> Cp + S
    r4 = CpS * k_d2
    # AtoC binds AtoS: C + S -> CS
    r5 = C * S * k_b3
    # AtoC unbinds AtoS: CS -> C + S
    r6 = CS * k_d3
    # Acetoacetate phosphorylates AtoS: S -> Sp
    r7 = S * k_ap
    # Dephosphorylation of AtoSP: Sp -> S
    r8 = Sp * k_ad
    # Phosphorylation of AtoC: CSp -> CpS
    r9 = CSp * k_pt
    # Dephosphorylation of AtoC: CpS -> CS (or should this be -> CSp?)
    r10 = CpS * k_ph
    # Binding of alternative dephosphatase: Cp + Ph -> PhCp
    r15 = Cp * Ph * k_b4
    # Unbinding of alternative dephosphatase: PhCp -> Cp + Ph
    r16 = PhCp * k_d4
    # Alternative dephosphorylation: PhCp -> C + Ph
    r17 = PhCp * k_cat
    # GFP degradation: GFP -> 0
    r20 = GFP * k_deg
    # Promoter binding: Cp + pato -> patoCp
    r21 = Cp * pato * k_bnd
    # Promoter unbinding: patoCp -> Cp + pato
    r22 = patoCp * k_unbnd
    # Promoter binding: C + pato -> patoC
    r23 = C * pato * k_nsbnd
    # Promoter unbinding: patoC -> C + pato
    r24 = patoC * k_nsunbnd
    # GFP transcription: patoCp -> patoCp + mGFP
    r25 = patoCp * k_pmgexp
    # GFP transcription: patoC -> patoC + mGFP
    r26 = patoC * k_npmgexp
    # Ribosome binding mRNA: mGFP + R -> mGFPR
    r27 = mGFP * R * k_mgbnd
    # GFP translation: mGFPR -> mGFP + R + uGFP
    r28 = mGFPR * k_gexp
    # GFP maturation: uGFP -> GFP
    r29 = uGFP * k_mat
    # mRNA degradation: mGFP -> 0
    r30 = mGFP * k_mgdeg

    ## ODEs
    dC = r2 + r6 + r17 + r24 - r1 - r5 - r23
    dS = r4 + r6 + r8 - r3 - r5 - r7
    dCp = r4 + r16 + r22 - r3 - r15 - r21
    dSp = r2 + r7 - r1 - r8
    dpato = r22 + r24 - r21 - r23
    dPh = r16 + r17 - r15
    dCS = r5 + r10 - r6
    dCSp = r1 - r2 - r9
    dCpS = r3 + r9 - r4 - r10
    dPhCp = r15 - r16 - r17
    dGFP = r29 - r20
    dpatoCp = r21 - r22
    dpatoC = r23 - r24
    dmGFP = r25 + r26 + r28 - r27 - r30
    dR = r28 - r27
    dmGFPR = r27 - r28
    duGFP = r28 - r29

    ds = (dC, dS, dCp, dSp, dpato, dPh, dCS, dCSp, dCpS, dPhCp, dpatoCp, dpatoC, dGFP, dmGFP, dR, dmGFPR, duGFP)

    return ds


def AtoSC_ori_d_dot(s, t, params):
    ## PARAMETERS
    k_b1, k_d1, k_b2, k_d2, k_b3, k_d3, k_ap, k_ad, k_pt, k_ph, k_dim, k_mon, k_b4, k_d4, k_cat, k_deg, k_dbnd, k_dunbnd, k_dmgexp, k_mgbnd, k_gexp, k_mat, k_mgdeg, k_lgexp = params

    ## SPECIES
    C, S, Cp, Sp, pato, Ph, Cpd, CS, CSp, CpS, patoCpd, PhCp, GFP, mGFP, R, mGFPR, uGFP = s

    ## REACTIONS
    # AtoC binds AtoSP: C + Sp -> CSp
    r1 = C * Sp * k_b1
    # AtoC unbinds AtoSP: CSp -> C + Sp
    r2 = CSp * k_d1
    # AtoCP binds AtoS: Cp + S -> CpS
    r3 = S * Cp * k_b2
    # AtoCP unbinds AtoS: CpS -> Cp + S
    r4 = CpS * k_d2
    # AtoC binds AtoS: C + S -> CS
    r5 = C * S * k_b3
    # AtoC unbinds AtoS: CS -> C + S
    r6 = CS * k_d3
    # Acetoacetate phosphorylates AtoS: S -> Sp
    r7 = S * k_ap
    # Dephosphorylation of AtoSP: Sp -> S
    r8 = Sp * k_ad
    # Phosphorylation of AtoC: CSp -> CpS
    r9 = CSp * k_pt
    # Dephosphorylation of AtoC: CpS -> CS (or should this be -> CSp?)
    r10 = CpS * k_ph
    # Dimerisation of AtoCP: Cp + Cp -> Cpd
    r11 = Cp * Cp * k_dim
    # Monomeristaion of AtoCP: Cpd -> Cp + Cp
    r12 = Cpd * k_mon
    # Promoter binding: pato + Cpd -> patoCpd
    r13 = pato * Cpd * k_dbnd
    # Promoter unbinding: patoCpd -> pato + Cpd
    r14 = patoCpd * k_dunbnd
    # Binding of alternative dephosphatase: Cp + Ph -> PhCp
    r15 = Cp * Ph * k_b4
    # Unbinding of alternative dephosphatase: PhCp -> Cp + Ph
    r16 = PhCp * k_d4
    # Alternative dephosphorylation: PhCp -> C + Ph
    r17 = PhCp * k_cat
    # GFP degradation: GFP -> 0
    r20 = GFP * k_deg
    # # Promoter binding: Cp + pato -> patoCp
    # r21 = Cp * pato * k_bnd
    # # Promoter unbinding: patoCp -> Cp + pato
    # r22 = patoCp * k_unbnd
    # # Promoter binding: C + pato -> patoC
    # r23 = C * pato * k_nsbnd
    # # Promoter unbinding: patoC -> C + pato
    # r24 = patoC * k_nsunbnd
    # # GFP transcription: patoCp -> patoCp + mGFP
    # r25 = patoCp * k_pmgexp
    # # GFP transcription: patoC -> patoC + mGFP
    # r26 = patoC * k_npmgexp
    # Ribosome binding mRNA: mGFP + R -> mGFPR
    r27 = mGFP * R * k_mgbnd
    # GFP translation: mGFPR -> mGFP + R + uGFP
    r28 = mGFPR * k_gexp
    # GFP maturation: uGFP -> GFP
    r29 = uGFP * k_mat
    # mRNA degradation: mGFP -> 0
    r30 = mGFP * k_mgdeg
    # GFP transcription: patoCpd -> patoCpd + mGFP
    r31 = patoCpd * k_dmgexp
    # leaky GFP transcription: pato -> pato + mGFP
    r32 = pato * k_lgexp

    ## ODEs
    dC = r2 + r6 + r17 - r1 - r5
    dS = r4 + r6 + r8 - r3 - r5 - r7
    dCp = r4 + 2 * r12 + r16 - r3 - 2 * r12 - r15
    dSp = r2 + r7 - r1 - r8
    dpato = r14 - r13
    dPh = r16 + r17 - r15
    dCpd = r11 + r14 - r12 - r13
    dCS = r5 + r10 - r6
    dCSp = r1 - r2 - r9
    dCpS = r3 + r9 - r4 - r10
    dpatoCpd = r13 - r14
    dPhCp = r15 - r16 - r17
    dGFP = r29 - r20
    dmGFP = r31 + r32 + r28 - r27 - r30
    dR = r28 - r27
    dmGFPR = r27 - r28
    duGFP = r28 - r29

    ds = (dC, dS, dCp, dSp, dpato, dPh, dCpd, dCS, dCSp, dCpS, dpatoCpd, dPhCp, dGFP, dmGFP, dR, dmGFPR, duGFP)

    return ds


def AtoSC_ori_s_dot(s, t, params):
    ## PARAMETERS
    k_b1, k_d1, k_b2, k_d2, k_b3, k_d3, k_ap, k_ad, k_pt, k_ph, k_deg, k_bnd, k_unbnd, k_pmgexp, k_gexp, k_mgdeg, k_lgexp = params

    ## SPECIES
    C, S, Cp, Sp, pato, CS, CSp, CpS, patoCp, GFP, mGFP = s

    ## REACTIONS
    # AtoC binds AtoSP: C + Sp -> CSp
    r1 = C * Sp * k_b1
    # AtoC unbinds AtoSP: CSp -> C + Sp
    r2 = CSp * k_d1
    # AtoCP binds AtoS: Cp + S -> CpS
    r3 = S * Cp * k_b2
    # AtoCP unbinds AtoS: CpS -> Cp + S
    r4 = CpS * k_d2
    # AtoC binds AtoS: C + S -> CS
    r5 = C * S * k_b3
    # AtoC unbinds AtoS: CS -> C + S
    r6 = CS * k_d3
    # Acetoacetate phosphorylates AtoS: S -> Sp
    r7 = S * k_ap
    # Dephosphorylation of AtoSP: Sp -> S
    r8 = Sp * k_ad
    # Phosphorylation of AtoC: CSp -> CpS
    r9 = CSp * k_pt
    # Dephosphorylation of AtoC: CpS -> CS (or should this be -> CSp?)
    r10 = CpS * k_ph
    # # Binding of alternative dephosphatase: Cp + Ph -> PhCp
    # r15 = Cp * Ph * k_b4
    # # Unbinding of alternative dephosphatase: PhCp -> Cp + Ph
    # r16 = PhCp * k_d4
    # # Alternative dephosphorylation: PhCp -> C + Ph
    # r17 = PhCp * k_cat
    # GFP degradation: GFP -> 0
    r20 = GFP * k_deg
    # Promoter binding: Cp + pato -> patoCp
    r21 = Cp * pato * k_bnd
    # Promoter unbinding: patoCp -> Cp + pato
    r22 = patoCp * k_unbnd
    # # Promoter binding: C + pato -> patoC
    # r23 = C * pato * k_nsbnd
    # # Promoter unbinding: patoC -> C + pato
    # r24 = patoC * k_nsunbnd
    # GFP transcription: patoCp -> patoCp + mGFP
    r25 = patoCp * k_pmgexp
    # # GFP transcription: patoC -> patoC + mGFP
    # r26 = patoC * k_npmgexp
    # # Ribosome binding mRNA: mGFP + R -> mGFPR
    # r27 = mGFP * R * k_mgbnd
    # # GFP translation: mGFPR -> mGFP + R + uGFP
    # r28 = mGFPR * k_gexp
    # # GFP maturation: uGFP -> GFP
    # r29 = uGFP * k_mat
    # mRNA degradation: mGFP -> 0
    r30 = mGFP * k_mgdeg
    # leaky GFP transcription: pato -> pato + mGFP
    r32 = pato * k_lgexp
    # # uGFP degredation: uGFP -> 0
    # r33 = uGFP * k_deg
    # GFP translation: mGFP -> mGFP + GFP
    r34 = mGFP * k_gexp

    ## ODEs
    dC = r2 + r6 - r1 - r5
    dS = r4 + r6 + r8 - r3 - r5 - r7
    dCp = r4 + r22 - r3 - r21
    dSp = r2 + r7 - r1 - r8
    dpato = r22 - r21
    dCS = r5 + r10 - r6
    dCSp = r1 - r2 - r9
    dCpS = r3 + r9 - r4 - r10
    dGFP = r34 - r20
    dpatoCp = r21 - r22
    dmGFP = r25 + r32 - r30

    ds = (dC, dS, dCp, dSp, dpato, dCS, dCSp, dCpS, dpatoCp, dGFP, dmGFP)

    return ds

