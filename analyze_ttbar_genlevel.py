#!/usr/bin/env python3
import os
import math
import gzip
import ROOT

ROOT.gROOT.SetBatch(True)


# =========================
# User paths (defaults)
# =========================
DEFAULT_LHEC = "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHeC_ttbar_semileptonic_decay/Events/run_01/LHeC_ttbar_semileptonic_decay.lhe"
DEFAULT_LHMUC = "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHmuC_ttbar_semileptonic_decay/Events/run_01/LHmuC_ttbar_semileptonic_decay.lhe"


# =========================
# Helpers
# =========================
def open_any(path):
    """Open plain text or .gz transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")


def safe_eta(px, py, pz):
    """Compute pseudorapidity from 3-momentum; robust near |p|=|pz|."""
    p2 = px*px + py*py + pz*pz
    p = math.sqrt(p2)
    # avoid division by zero
    if p == abs(pz):
        return 1e9 if pz >= 0 else -1e9
    return 0.5 * math.log((p + pz) / (p - pz))


def pt(px, py):
    return math.sqrt(px*px + py*py)


class Particle:
    __slots__ = ("idx", "pid", "st", "m1", "m2", "px", "py", "pz", "E", "m")
    def __init__(self, idx, pid, st, m1, m2, px, py, pz, E, m):
        self.idx = idx          # 1-based LHE index inside event
        self.pid = pid
        self.st  = st
        self.m1  = m1
        self.m2  = m2
        self.px  = px
        self.py  = py
        self.pz  = pz
        self.E   = E
        self.m   = m

    def eta(self):
        return safe_eta(self.px, self.py, self.pz)

    def pt(self):
        return pt(self.px, self.py)

    def p4(self):
        v = ROOT.TLorentzVector()
        v.SetPxPyPzE(self.px, self.py, self.pz, self.E)
        return v


def children(parts, parent_idx):
    return [p for p in parts if p.m1 == parent_idx or p.m2 == parent_idx]


# =========================
# Read xsec and init info
# =========================
def read_mg_generation_info(lhe_path):
    """
    Parse <MGGenerationInfo> block:
      Number of Events
      Integrated weight (pb)
    Returns: (nevents, xsec_pb) or (None, None) if not found.
    """
    nevt = None
    xsec = None
    inside = False
    with open_any(lhe_path) as f:
        for line in f:
            if "<MGGenerationInfo>" in line:
                inside = True
                continue
            if inside and "</MGGenerationInfo>" in line:
                break
            if inside:
                if "Number of Events" in line:
                    # ... : 1000000
                    nevt = int(line.split(":")[-1].strip())
                if "Integrated weight (pb)" in line:
                    xsec = float(line.split(":")[-1].strip())
    return nevt, xsec


def read_init_subprocess_xsecs(lhe_path):
    """
    Parse <init> block and return:
      - beams (pdg1,pdg2,E1,E2)
      - dict proc_index -> xsec_pb (XSECUP)
      - total_xsec_pb = sum(XSECUP)
    """
    inside = False
    got_first = False
    nprup = None
    beams = None
    xsecs = {}
    with open_any(lhe_path) as f:
        for line in f:
            if "<init>" in line:
                inside = True
                continue
            if inside and "</init>" in line:
                break
            if not inside:
                continue
            s = line.strip()
            if not s:
                continue
            if not got_first:
                # first init line: ... NPRUP at end
                parts = s.split()
                # beam PDGs and energies
                pdg1 = int(parts[0]); pdg2 = int(parts[1])
                E1 = float(parts[2]); E2 = float(parts[3])
                nprup = int(parts[-1])
                beams = (pdg1, pdg2, E1, E2)
                got_first = True
                continue
            # next nprup lines: XSECUP XERRUP XMAXUP LPRUP
            if nprup is not None and len(xsecs) < nprup:
                cols = s.split()
                xsecup = float(cols[0])
                lprup  = int(cols[3])
                # MG typically uses LPRUP as subprocess label; but IDPRUP is 1..NPRUP
                # We'll store in read order (1..NPRUP)
                idx = len(xsecs) + 1
                xsecs[idx] = xsecup
    total = sum(xsecs.values()) if xsecs else None
    return beams, xsecs, total


# =========================
# Stream events from LHE
# =========================
def iter_event_blocks(lhe_path):
    """Yield list of lines (stripped) inside each <event>...</event>."""
    with open_any(lhe_path) as f:
        inside = False
        buf = []
        for line in f:
            if "<event>" in line:
                inside = True
                buf = []
                continue
            if "</event>" in line and inside:
                yield buf
                inside = False
                buf = []
                continue
            if inside:
                s = line.strip()
                if s:
                    buf.append(s)


def parse_event(block_lines):
    """
    block_lines[0] = event header: NUP IDPRUP XWGTUP SCALUP AQEDUP AQCDUP
    next NUP lines are particles.
    returns: (idprup, weight, parts_list)
    """
    header = block_lines[0].split()
    nup = int(header[0])
    idprup = int(header[1])
    wgt = float(header[2])

    parts = []
    for i in range(1, 1 + nup):
        cols = block_lines[i].split()
        pid = int(cols[0])
        st  = int(cols[1])
        m1  = int(cols[2])
        m2  = int(cols[3])
        px  = float(cols[6])
        py  = float(cols[7])
        pz  = float(cols[8])
        E   = float(cols[9])
        m   = float(cols[10])
        parts.append(Particle(idx=i, pid=pid, st=st, m1=m1, m2=m2, px=px, py=py, pz=pz, E=E, m=m))
    return idprup, wgt, parts


# =========================
# Reconstruct ttbar semileptonic
# =========================
LEP_IDS = {11, -11, 13, -13, 15, -15}
NU_IDS  = {12, -12, 14, -14, 16, -16}
QUARK_IDS = {1, -1, 2, -2, 3, -3, 4, -4, 5, -5}

def reconstruct_ttbar(parts):
    """
    Reconstruct:
      - decay lepton (from W)
      - 4 decay quarks (b, bbar, q, q')
      - leptonic top 4-vector
      - hadronic top 4-vector
      - hadronic W 4-vector
    Uses mother pointers (m1/m2).
    Returns dict or None if can't reconstruct.
    """

    # Find top and antitop (status 2 is typical in MG LHE)
    t    = next((p for p in parts if p.pid == 6), None)
    tbar = next((p for p in parts if p.pid == -6), None)
    if t is None or tbar is None:
        return None

    def get_top_decay(top_p):
        ch = children(parts, top_p.idx)
        # identify b and W among children
        b = next((x for x in ch if x.pid == (5 if top_p.pid == 6 else -5)), None)
        W = next((x for x in ch if x.pid == (24 if top_p.pid == 6 else -24)), None)
        # fallback if exact sign not found (rare)
        if b is None:
            b = next((x for x in ch if abs(x.pid) == 5), None)
        if W is None:
            W = next((x for x in ch if abs(x.pid) == 24), None)
        if b is None or W is None:
            return None
        wch = children(parts, W.idx)
        return b, W, wch

    dec_t = get_top_decay(t)
    dec_tb = get_top_decay(tbar)
    if dec_t is None or dec_tb is None:
        return None

    b_t,  Wp,  Wp_ch  = dec_t
    b_tb, Wm,  Wm_ch  = dec_tb

    def is_leptonic(wch):
        has_lep = any(abs(x.pid) in {11,13,15} for x in wch)
        has_nu  = any(abs(x.pid) in {12,14,16} for x in wch)
        return has_lep and has_nu

    # Determine which W is leptonic
    if is_leptonic(Wp_ch):
        lepW_ch = Wp_ch
        hadW_ch = Wm_ch
        b_lep = b_t
        b_had = b_tb
    elif is_leptonic(Wm_ch):
        lepW_ch = Wm_ch
        hadW_ch = Wp_ch
        b_lep = b_tb
        b_had = b_t
    else:
        return None

    lepton = next((x for x in lepW_ch if x.pid in LEP_IDS and x.st == 1), None)
    nu     = next((x for x in lepW_ch if x.pid in NU_IDS and x.st == 1), None)
    if lepton is None or nu is None:
        return None

    # Hadronic W: pick quarks (status 1)
    qq = [x for x in hadW_ch if x.pid in QUARK_IDS and x.st == 1]
    if len(qq) != 2:
        # Sometimes W children might include something unexpected; fail safely
        return None

    q1, q2 = qq[0], qq[1]

    # Reconstruct 4-vectors using TLorentzVector
    p4_top_lep = b_lep.p4() + lepton.p4() + nu.p4()
    p4_top_had = b_had.p4() + q1.p4() + q2.p4()
    p4_W_had   = q1.p4() + q2.p4()

    jets = [b_lep, b_had, q1, q2]

    return {
        "lepton": lepton,
        "nu": nu,
        "jets": jets,
        "bjets": [b_lep, b_had],
        "ljets": [q1, q2],
        "top_lep": p4_top_lep,
        "top_had": p4_top_had,
        "w_had": p4_W_had,
    }


# =========================
# Histogram booking
# =========================
def book_hists(tag):
    """
    Create two sets: inclusive and fiducial (after eta cuts).
    """
    h = {}
    def H1(name, title, nb, x1, x2):
        hh = ROOT.TH1F(f"{tag}_{name}", f"{title};{title.split(';')[1] if ';' in title else ''}", nb, x1, x2)
        hh.Sumw2()
        return hh

    # Inclusive
    h["lep_eta"]   = ROOT.TH1F(f"{tag}_lep_eta",  "lepton eta;eta;events", 120, -6, 6)
    h["lep_pt"]    = ROOT.TH1F(f"{tag}_lep_pt",   "lepton pT;pT [GeV];events", 200, 0, 1000)

    h["jet_eta"]   = ROOT.TH1F(f"{tag}_jet_eta",  "decay-jet eta;eta;jets", 120, -6, 6)
    h["jet_pt"]    = ROOT.TH1F(f"{tag}_jet_pt",   "decay-jet pT;pT [GeV];jets", 200, 0, 1000)

    h["b_eta"]     = ROOT.TH1F(f"{tag}_b_eta",    "b-jet eta;eta;b-jets", 120, -6, 6)
    h["b_pt"]      = ROOT.TH1F(f"{tag}_b_pt",     "b-jet pT;pT [GeV];b-jets", 200, 0, 1000)

    h["lj_eta"]    = ROOT.TH1F(f"{tag}_lj_eta",   "light-jet eta;eta;light-jets", 120, -6, 6)
    h["lj_pt"]     = ROOT.TH1F(f"{tag}_lj_pt",    "light-jet pT;pT [GeV];light-jets", 200, 0, 1000)

    h["topL_m"]    = ROOT.TH1F(f"{tag}_topL_m",   "m(top_lep);m [GeV];events", 200, 0, 400)
    h["topH_m"]    = ROOT.TH1F(f"{tag}_topH_m",   "m(top_had);m [GeV];events", 200, 0, 400)
    h["topL_pt"]   = ROOT.TH1F(f"{tag}_topL_pt",  "pT(top_lep);pT [GeV];events", 200, 0, 2000)
    h["topH_pt"]   = ROOT.TH1F(f"{tag}_topH_pt",  "pT(top_had);pT [GeV];events", 200, 0, 2000)
    h["topL_eta"]  = ROOT.TH1F(f"{tag}_topL_eta", "eta(top_lep);eta;events", 120, -6, 6)
    h["topH_eta"]  = ROOT.TH1F(f"{tag}_topH_eta", "eta(top_had);eta;events", 120, -6, 6)

    h["Whad_m"]    = ROOT.TH1F(f"{tag}_Whad_m",   "m(W_had);m [GeV];events", 160, 0, 200)

    # Fiducial duplicates
    hf = {}
    for k, hh in h.items():
        name = hh.GetName().replace(tag, tag + "_fid")
        title = hh.GetTitle() + " (fiducial)"
        nb = hh.GetNbinsX()
        x1 = hh.GetXaxis().GetXmin()
        x2 = hh.GetXaxis().GetXmax()
        hnew = ROOT.TH1F(name, title + f";{hh.GetXaxis().GetTitle()};{hh.GetYaxis().GetTitle()}", nb, x1, x2)
        hnew.Sumw2()
        hf[k] = hnew

    return h, hf


def save_pngs(rootfile, outdir):
    """
    Save all TH1 in a ROOT file as PNGs.
    """
    os.makedirs(outdir, exist_ok=True)
    f = ROOT.TFile.Open(rootfile, "READ")
    keys = f.GetListOfKeys()
    c = ROOT.TCanvas("c", "c", 900, 700)
    for k in keys:
        obj = k.ReadObj()
        if not obj.InheritsFrom("TH1"):
            continue
        obj.SetLineWidth(2)
        obj.Draw("HIST")
        name = obj.GetName()
        c.SaveAs(os.path.join(outdir, f"{name}.png"))
    f.Close()


# =========================
# Core analysis driver
# =========================
def analyze_file(lhe_path, eta_max, tag, outdir, max_events=None):
    """
    Analyze one LHE file and produce:
      - ROOT hist file
      - PNGs
      - summary text
    Returns dict summary.
    """
    os.makedirs(outdir, exist_ok=True)

    nevt_mg, xsec_mg = read_mg_generation_info(lhe_path)
    beams, xsecs_proc, xsec_init_sum = read_init_subprocess_xsecs(lhe_path)

    # Use MGGenerationInfo xsec if available, otherwise init sum
    xsec_total = xsec_mg if xsec_mg is not None else xsec_init_sum

    # Book histograms
    h_all, h_fid = book_hists(tag)

    # Counters per subprocess IDPRUP
    # weighted sums (using XWGTUP)
    sumw_all = {}
    sumw_pass = {}
    nev_all = {}
    nev_pass = {}

    # Total counters
    sumw_all_tot = 0.0
    sumw_pass_tot = 0.0
    nev_all_tot = 0
    nev_pass_tot = 0

    nev_read = 0
    nev_reco = 0

    for block in iter_event_blocks(lhe_path):
        idprup, wgt, parts = parse_event(block)
        nev_read += 1
        if max_events and nev_read > max_events:
            break

        rec = reconstruct_ttbar(parts)
        if rec is None:
            continue
        nev_reco += 1

        # init dicts
        sumw_all.setdefault(idprup, 0.0)
        sumw_pass.setdefault(idprup, 0.0)
        nev_all.setdefault(idprup, 0)
        nev_pass.setdefault(idprup, 0)

        sumw_all[idprup] += wgt
        nev_all[idprup] += 1
        sumw_all_tot += wgt
        nev_all_tot += 1

        lep = rec["lepton"]
        jets = rec["jets"]

        # Fill inclusive hists (weight by wgt)
        h_all["lep_eta"].Fill(lep.eta(), wgt)
        h_all["lep_pt"].Fill(lep.pt(), wgt)

        # sort jets by pT for some “leading jet” behavior (still fill all)
        jets_sorted = sorted(jets, key=lambda x: x.pt(), reverse=True)

        for j in jets:
            h_all["jet_eta"].Fill(j.eta(), wgt)
            h_all["jet_pt"].Fill(j.pt(), wgt)

        for b in rec["bjets"]:
            h_all["b_eta"].Fill(b.eta(), wgt)
            h_all["b_pt"].Fill(b.pt(), wgt)

        for lj in rec["ljets"]:
            h_all["lj_eta"].Fill(lj.eta(), wgt)
            h_all["lj_pt"].Fill(lj.pt(), wgt)

        # tops
        h_all["topL_m"].Fill(rec["top_lep"].M(), wgt)
        h_all["topH_m"].Fill(rec["top_had"].M(), wgt)
        h_all["topL_pt"].Fill(rec["top_lep"].Pt(), wgt)
        h_all["topH_pt"].Fill(rec["top_had"].Pt(), wgt)
        h_all["topL_eta"].Fill(rec["top_lep"].Eta(), wgt)
        h_all["topH_eta"].Fill(rec["top_had"].Eta(), wgt)
        h_all["Whad_m"].Fill(rec["w_had"].M(), wgt)

        # Krzysztof GEN-level fiducial:
        pass_eta = (abs(lep.eta()) < eta_max) and all(abs(j.eta()) < eta_max for j in jets)

        if pass_eta:
            sumw_pass[idprup] += wgt
            nev_pass[idprup] += 1
            sumw_pass_tot += wgt
            nev_pass_tot += 1

            # Fill fiducial hists
            h_fid["lep_eta"].Fill(lep.eta(), wgt)
            h_fid["lep_pt"].Fill(lep.pt(), wgt)

            for j in jets:
                h_fid["jet_eta"].Fill(j.eta(), wgt)
                h_fid["jet_pt"].Fill(j.pt(), wgt)

            for b in rec["bjets"]:
                h_fid["b_eta"].Fill(b.eta(), wgt)
                h_fid["b_pt"].Fill(b.pt(), wgt)

            for lj in rec["ljets"]:
                h_fid["lj_eta"].Fill(lj.eta(), wgt)
                h_fid["lj_pt"].Fill(lj.pt(), wgt)

            h_fid["topL_m"].Fill(rec["top_lep"].M(), wgt)
            h_fid["topH_m"].Fill(rec["top_had"].M(), wgt)
            h_fid["topL_pt"].Fill(rec["top_lep"].Pt(), wgt)
            h_fid["topH_pt"].Fill(rec["top_had"].Pt(), wgt)
            h_fid["topL_eta"].Fill(rec["top_lep"].Eta(), wgt)
            h_fid["topH_eta"].Fill(rec["top_had"].Eta(), wgt)
            h_fid["Whad_m"].Fill(rec["w_had"].M(), wgt)

    # -------------------------
    # Compute fiducial xsec
    # Best practice: per-subprocess
    # sigma_fid = sum_i XSECUP_i * A_i  where A_i = sumw_pass_i/sumw_all_i
    # -------------------------
    sigma_fid = 0.0
    perproc_lines = []
    for iproc, xsec_i in xsecs_proc.items():
        sw_all = sumw_all.get(iproc, 0.0)
        sw_pas = sumw_pass.get(iproc, 0.0)
        Ai = (sw_pas / sw_all) if sw_all > 0 else 0.0
        sigma_i_fid = xsec_i * Ai
        sigma_fid += sigma_i_fid

        perproc_lines.append(
            f"  proc {iproc:2d}: XSECUP={xsec_i:.12g} pb | sumw_all={sw_all:.6g} | sumw_pass={sw_pas:.6g} | A={Ai:.6g} | sigma_fid_i={sigma_i_fid:.12g} pb"
        )

    # Also compute total-weight acceptance (cross-check)
    A_tot = (sumw_pass_tot / sumw_all_tot) if sumw_all_tot > 0 else 0.0
    sigma_fid_alt = (xsec_total * A_tot) if xsec_total is not None else None

    # -------------------------
    # Write ROOT file
    # -------------------------
    root_out = os.path.join(outdir, f"{tag}_hist.root")
    fout = ROOT.TFile(root_out, "RECREATE")
    for hh in h_all.values():
        hh.Write()
    for hh in h_fid.values():
        hh.Write()
    fout.Close()

    # Save PNGs
    png_dir = os.path.join(outdir, "png")
    save_pngs(root_out, png_dir)

    # Summary text
    txt_out = os.path.join(outdir, f"{tag}_summary.txt")
    with open(txt_out, "w") as fp:
        fp.write(f"File: {lhe_path}\n")
        fp.write(f"Tag: {tag}\n")
        fp.write(f"Beams (pdg1,pdg2,E1,E2): {beams}\n")
        fp.write(f"MGGenerationInfo: nevents={nevt_mg}, xsec={xsec_mg} pb\n")
        fp.write(f"Init xsec sum: {xsec_init_sum} pb\n")
        fp.write(f"Using total xsec: {xsec_total} pb\n")
        fp.write(f"Events read: {nev_read}\n")
        fp.write(f"Events reconstructed (ttbar ok): {nev_reco}\n")
        fp.write(f"Fiducial cut: |eta_lep|<{eta_max} AND all 4 decay jets |eta|<{eta_max}\n")
        fp.write(f"Total weighted acceptance A_tot = {A_tot:.12g}\n")
        fp.write(f"Sigma_fid (per-proc) = {sigma_fid:.12g} pb\n")
        if sigma_fid_alt is not None:
            fp.write(f"Sigma_fid_alt (total*Atot) = {sigma_fid_alt:.12g} pb\n")
        fp.write("\nPer-subprocess:\n")
        fp.write("\n".join(perproc_lines) + "\n")

    # Print short summary
    print(f"\n=== {tag} ===")
    print(f"File: {lhe_path}")
    print(f"Beams: {beams}")
    print(f"Total xsec: {xsec_total} pb")
    print(f"Events read: {nev_read}, reconstructed: {nev_reco}")
    print(f"Fiducial: |eta|<{eta_max} on lepton + all 4 decay jets")
    print(f"A_tot = {A_tot:.6g}")
    print(f"sigma_fid (per-proc) = {sigma_fid:.12g} pb")
    if sigma_fid_alt is not None:
        print(f"sigma_fid_alt (total*Atot) = {sigma_fid_alt:.12g} pb")
    for line in perproc_lines:
        print(line)

    return {
        "tag": tag,
        "lhe": lhe_path,
        "beams": beams,
        "xsec_total_pb": xsec_total,
        "sigma_fid_pb": sigma_fid,
        "A_tot": A_tot,
        "root_out": root_out,
        "png_dir": png_dir,
        "txt_out": txt_out,
    }


# =========================
# Main: analyze both samples
# =========================
def main():
    import argparse
    ap = argparse.ArgumentParser(description="GEN-level ttbar semileptonic analysis for LHeC/LHmuC LHE files.")
    ap.add_argument("--lhec", default=DEFAULT_LHEC, help="Path to LHeC LHE file")
    ap.add_argument("--lhmuc", default=DEFAULT_LHMUC, help="Path to LHmuC LHE file")
    ap.add_argument("--outdir", default="ttbar_genlevel_outputs", help="Output directory")
    ap.add_argument("--maxev", type=int, default=None, help="Optional: max number of events to process")
    args = ap.parse_args()

    # LHeC: eta<5
    analyze_file(args.lhec, eta_max=5.0, tag="LHeC", outdir=os.path.join(args.outdir, "LHeC"), max_events=args.maxev)

    # LHmuC: eta<4
    analyze_file(args.lhmuc, eta_max=4.0, tag="LHmuC", outdir=os.path.join(args.outdir, "LHmuC"), max_events=args.maxev)


if __name__ == "__main__":
    main()

