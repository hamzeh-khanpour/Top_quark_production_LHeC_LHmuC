
#!/usr/bin/env python3
import os
import math
import gzip
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt


# -------------------------
# Your files (as requested)
# -------------------------
LHEC_LHE = "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHeC_ttbar_semileptonic_decay/Events/run_01/LHeC_ttbar_semileptonic_decay.lhe"
LHMUC_LHE = "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHmuC_ttbar_semileptonic_decay/Events/run_01/LHmuC_ttbar_semileptonic_decay.lhe"


def open_any(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")


# -------------------------
# Simple Lorentz vector
# -------------------------
class LV:
    __slots__ = ("px", "py", "pz", "E")
    def __init__(self, px, py, pz, E):
        self.px, self.py, self.pz, self.E = px, py, pz, E

    def __add__(self, other):
        return LV(self.px + other.px, self.py + other.py, self.pz + other.pz, self.E + other.E)

    def pt(self):
        return math.hypot(self.px, self.py)

    def p(self):
        return math.sqrt(self.px*self.px + self.py*self.py + self.pz*self.pz)

    def eta(self):
        p = self.p()
        if p == abs(self.pz):
            return 1e9 if self.pz >= 0 else -1e9
        return 0.5 * math.log((p + self.pz) / (p - self.pz))

    def m2(self):
        return self.E*self.E - (self.px*self.px + self.py*self.py + self.pz*self.pz)

    def m(self):
        m2 = self.m2()
        return math.sqrt(m2) if m2 > 0 else 0.0


# -------------------------
# Particle record
# -------------------------
class P:
    __slots__ = ("idx", "pid", "st", "m1", "m2", "lv")
    def __init__(self, idx, pid, st, m1, m2, px, py, pz, E):
        self.idx = idx  # 1-based
        self.pid = pid
        self.st = st
        self.m1 = m1
        self.m2 = m2
        self.lv = LV(px, py, pz, E)


def children(parts, parent_idx):
    # mother pointers are 1-based indices; 0 means none
    return [p for p in parts if p.m1 == parent_idx or p.m2 == parent_idx]


# -------------------------
# Read xsec blocks
# -------------------------
def read_mg_generation_info(lhe_path):
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
                    nevt = int(line.split(":")[-1].strip())
                elif "Integrated weight (pb)" in line:
                    xsec = float(line.split(":")[-1].strip())
    return nevt, xsec


def read_init_subprocess_xsecs(lhe_path):
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
                cols = s.split()
                pdg1, pdg2 = int(cols[0]), int(cols[1])
                E1, E2 = float(cols[2]), float(cols[3])
                nprup = int(cols[-1])
                beams = (pdg1, pdg2, E1, E2, nprup)
                got_first = True
                continue
            if nprup is not None and len(xsecs) < nprup:
                cols = s.split()
                xsecup = float(cols[0])
                iproc = len(xsecs) + 1  # by read order, matches IDPRUP=1..NPRUP
                xsecs[iproc] = xsecup
    return beams, xsecs, sum(xsecs.values()) if xsecs else None


# -------------------------
# Event streaming
# -------------------------
def iter_event_blocks(lhe_path):
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


def parse_event(block):
    # header: NUP IDPRUP XWGTUP SCALUP AQEDUP AQCDUP
    hdr = block[0].split()
    nup = int(hdr[0])
    idprup = int(hdr[1])
    wgt = float(hdr[2])

    parts = []
    for i in range(1, 1 + nup):
        cols = block[i].split()
        pid = int(cols[0])
        st = int(cols[1])
        m1 = int(cols[2])
        m2 = int(cols[3])
        px = float(cols[6])
        py = float(cols[7])
        pz = float(cols[8])
        E = float(cols[9])
        parts.append(P(i, pid, st, m1, m2, px, py, pz, E))

    return idprup, wgt, parts


# -------------------------
# Reconstruction (semileptonic ttbar)
# -------------------------
LEP = {11, -11, 13, -13, 15, -15}
NU  = {12, -12, 14, -14, 16, -16}
# for hadronic W we prefer light quarks (<=4) if possible
LIGHTQ = {1, -1, 2, -2, 3, -3, 4, -4}
QUARKS = {1, -1, 2, -2, 3, -3, 4, -4, 5, -5}

def reconstruct(parts):
    t = next((p for p in parts if p.pid == 6), None)
    tb = next((p for p in parts if p.pid == -6), None)
    if t is None or tb is None:
        return None

    def decay_of_top(top):
        ch = children(parts, top.idx)
        b = next((x for x in ch if abs(x.pid) == 5), None)
        W = next((x for x in ch if abs(x.pid) == 24), None)
        if b is None or W is None:
            return None
        wch = children(parts, W.idx)
        return b, W, wch

    dt = decay_of_top(t)
    dtb = decay_of_top(tb)
    if dt is None or dtb is None:
        return None

    b_t, W_t, Wt_ch = dt
    b_tb, W_tb, Wtb_ch = dtb

    def is_leptonic(wch):
        return any(x.pid in LEP and x.st == 1 for x in wch) and any(x.pid in NU and x.st == 1 for x in wch)

    if is_leptonic(Wt_ch):
        lep_ch, had_ch = Wt_ch, Wtb_ch
        b_lep, b_had = b_t, b_tb
    elif is_leptonic(Wtb_ch):
        lep_ch, had_ch = Wtb_ch, Wt_ch
        b_lep, b_had = b_tb, b_t
    else:
        return None

    lepton = next((x for x in lep_ch if x.pid in LEP and x.st == 1), None)
    nu = next((x for x in lep_ch if x.pid in NU and x.st == 1), None)
    if lepton is None or nu is None:
        return None

    # hadronic W quarks: prefer light (no b) if available
    q_light = [x for x in had_ch if x.pid in LIGHTQ and x.st == 1]
    if len(q_light) == 2:
        q1, q2 = q_light
    else:
        q_any = [x for x in had_ch if x.pid in QUARKS and x.st == 1]
        if len(q_any) != 2:
            return None
        q1, q2 = q_any

    jets = [b_lep, b_had, q1, q2]

    top_lep = b_lep.lv + lepton.lv + nu.lv
    top_had = b_had.lv + q1.lv + q2.lv
    w_had = q1.lv + q2.lv

    return {
        "id_lep": lepton,
        "jets": jets,
        "bjets": [b_lep, b_had],
        "ljets": [q1, q2],
        "top_lep": top_lep,
        "top_had": top_had,
        "w_had": w_had,
    }


# -------------------------
# Histogram helper (weighted)
# -------------------------
class Hist1D:
    def __init__(self, name, bins):
        self.name = name
        self.bins = np.asarray(bins, dtype=float)
        self.counts = np.zeros(len(self.bins) - 1, dtype=float)

    def fill(self, x, w=1.0):
        # under/overflow ignored (can be added if you want)
        i = np.searchsorted(self.bins, x, side="right") - 1
        if 0 <= i < len(self.counts):
            self.counts[i] += w

    def save_png(self, outpath, xlabel):
        centers = 0.5 * (self.bins[:-1] + self.bins[1:])
        widths = self.bins[1:] - self.bins[:-1]
        plt.figure()
        plt.bar(centers, self.counts, width=widths, align="center")
        plt.xlabel(xlabel)
        plt.ylabel("weighted events")
        plt.title(self.name)
        plt.tight_layout()
        plt.savefig(outpath)
        plt.close()


def book_hists(prefix):
    # binning similar to what we discussed
    bins_eta = np.linspace(-6, 6, 121)
    bins_pt1 = np.linspace(0, 1000, 201)
    bins_pt2 = np.linspace(0, 2000, 201)
    bins_mtop = np.linspace(0, 400, 201)
    bins_mw = np.linspace(0, 200, 161)

    h = {}
    def add(n, b): h[n] = Hist1D(f"{prefix}_{n}", b)

    add("lep_eta", bins_eta)
    add("lep_pt", bins_pt1)
    add("jet_eta", bins_eta)
    add("jet_pt", bins_pt1)
    add("b_eta", bins_eta)
    add("b_pt", bins_pt1)
    add("lj_eta", bins_eta)
    add("lj_pt", bins_pt1)
    add("topL_eta", bins_eta)
    add("topH_eta", bins_eta)
    add("topL_pt", bins_pt2)
    add("topH_pt", bins_pt2)
    add("topL_m", bins_mtop)
    add("topH_m", bins_mtop)
    add("Whad_m", bins_mw)

    return h


# -------------------------
# Analysis for one sample
# -------------------------
def analyze_sample(lhe_path, outdir, tag, max_events=None):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "png"), exist_ok=True)

    nevt_mg, xsec_mg = read_mg_generation_info(lhe_path)
    beams, xsecs_proc, xsec_init = read_init_subprocess_xsecs(lhe_path)

    xsec_total = xsec_mg if xsec_mg is not None else xsec_init

    # Histograms: inclusive + fid(eta<5) + fid(eta<4)
    h_all = book_hists(f"{tag}_ALL")
    h_f5  = book_hists(f"{tag}_FID5")
    h_f4  = book_hists(f"{tag}_FID4")

    # weighted acceptance counters per subprocess for eta<5 and eta<4
    sw_all = defaultdict(float)
    sw_p5  = defaultdict(float)
    sw_p4  = defaultdict(float)

    nev_read = 0
    nev_reco = 0

    for block in iter_event_blocks(lhe_path):
        nev_read += 1
        if max_events and nev_read > max_events:
            break

        idprup, wgt, parts = parse_event(block)
        rec = reconstruct(parts)
        if rec is None:
            continue
        nev_reco += 1

        sw_all[idprup] += wgt

        lep = rec["id_lep"]
        jets = rec["jets"]

        # inclusive fills
        h_all["lep_eta"].fill(lep.lv.eta(), wgt)
        h_all["lep_pt"].fill(lep.lv.pt(), wgt)
        for j in jets:
            h_all["jet_eta"].fill(j.lv.eta(), wgt)
            h_all["jet_pt"].fill(j.lv.pt(), wgt)
        for b in rec["bjets"]:
            h_all["b_eta"].fill(b.lv.eta(), wgt)
            h_all["b_pt"].fill(b.lv.pt(), wgt)
        for lj in rec["ljets"]:
            h_all["lj_eta"].fill(lj.lv.eta(), wgt)
            h_all["lj_pt"].fill(lj.lv.pt(), wgt)

        h_all["topL_eta"].fill(rec["top_lep"].eta(), wgt)
        h_all["topH_eta"].fill(rec["top_had"].eta(), wgt)
        h_all["topL_pt"].fill(rec["top_lep"].pt(), wgt)
        h_all["topH_pt"].fill(rec["top_had"].pt(), wgt)
        h_all["topL_m"].fill(rec["top_lep"].m(), wgt)
        h_all["topH_m"].fill(rec["top_had"].m(), wgt)
        h_all["Whad_m"].fill(rec["w_had"].m(), wgt)

        # fiducial cuts
        pass5 = (abs(lep.lv.eta()) < 5.0) and all(abs(j.lv.eta()) < 5.0 for j in jets)
        pass4 = (abs(lep.lv.eta()) < 4.0) and all(abs(j.lv.eta()) < 4.0 for j in jets)

        if pass5:
            sw_p5[idprup] += wgt
            # fill fid5 hists
            h_f5["lep_eta"].fill(lep.lv.eta(), wgt)
            h_f5["lep_pt"].fill(lep.lv.pt(), wgt)
            for j in jets:
                h_f5["jet_eta"].fill(j.lv.eta(), wgt)
                h_f5["jet_pt"].fill(j.lv.pt(), wgt)
            for b in rec["bjets"]:
                h_f5["b_eta"].fill(b.lv.eta(), wgt)
                h_f5["b_pt"].fill(b.lv.pt(), wgt)
            for lj in rec["ljets"]:
                h_f5["lj_eta"].fill(lj.lv.eta(), wgt)
                h_f5["lj_pt"].fill(lj.lv.pt(), wgt)
            h_f5["topL_eta"].fill(rec["top_lep"].eta(), wgt)
            h_f5["topH_eta"].fill(rec["top_had"].eta(), wgt)
            h_f5["topL_pt"].fill(rec["top_lep"].pt(), wgt)
            h_f5["topH_pt"].fill(rec["top_had"].pt(), wgt)
            h_f5["topL_m"].fill(rec["top_lep"].m(), wgt)
            h_f5["topH_m"].fill(rec["top_had"].m(), wgt)
            h_f5["Whad_m"].fill(rec["w_had"].m(), wgt)

        if pass4:
            sw_p4[idprup] += wgt
            # fill fid4 hists
            h_f4["lep_eta"].fill(lep.lv.eta(), wgt)
            h_f4["lep_pt"].fill(lep.lv.pt(), wgt)
            for j in jets:
                h_f4["jet_eta"].fill(j.lv.eta(), wgt)
                h_f4["jet_pt"].fill(j.lv.pt(), wgt)
            for b in rec["bjets"]:
                h_f4["b_eta"].fill(b.lv.eta(), wgt)
                h_f4["b_pt"].fill(b.lv.pt(), wgt)
            for lj in rec["ljets"]:
                h_f4["lj_eta"].fill(lj.lv.eta(), wgt)
                h_f4["lj_pt"].fill(lj.lv.pt(), wgt)
            h_f4["topL_eta"].fill(rec["top_lep"].eta(), wgt)
            h_f4["topH_eta"].fill(rec["top_had"].eta(), wgt)
            h_f4["topL_pt"].fill(rec["top_lep"].pt(), wgt)
            h_f4["topH_pt"].fill(rec["top_had"].pt(), wgt)
            h_f4["topL_m"].fill(rec["top_lep"].m(), wgt)
            h_f4["topH_m"].fill(rec["top_had"].m(), wgt)
            h_f4["Whad_m"].fill(rec["w_had"].m(), wgt)

    # Compute sigma_fid via per-subprocess formula:
    # sigma_fid = sum_i XSECUP_i * (sumw_pass_i / sumw_all_i)
    def sigma_fid(sw_pass):
        sig = 0.0
        for i, xsec_i in xsecs_proc.items():
            a = (sw_pass[i] / sw_all[i]) if sw_all[i] > 0 else 0.0
            sig += xsec_i * a
        return sig

    sigma_fid5 = sigma_fid(sw_p5)
    sigma_fid4 = sigma_fid(sw_p4)

    # Also total-weight acceptance for cross-check
    sw_all_tot = sum(sw_all.values())
    A5 = (sum(sw_p5.values()) / sw_all_tot) if sw_all_tot > 0 else 0.0
    A4 = (sum(sw_p4.values()) / sw_all_tot) if sw_all_tot > 0 else 0.0
    sigma_fid5_alt = xsec_total * A5 if xsec_total is not None else None
    sigma_fid4_alt = xsec_total * A4 if xsec_total is not None else None

    # Save histograms to npz
    npz_path = os.path.join(outdir, f"{tag}_hists.npz")
    def pack_histdict(hdict):
        out = {}
        for k, hh in hdict.items():
            out[f"{k}_bins"] = hh.bins
            out[f"{k}_counts"] = hh.counts
        return out

    np.savez_compressed(npz_path,
                        **{f"ALL_{k}": v for k, v in pack_histdict(h_all).items()},
                        **{f"FID5_{k}": v for k, v in pack_histdict(h_f5).items()},
                        **{f"FID4_{k}": v for k, v in pack_histdict(h_f4).items()},
                        )

    # Save PNGs
    def dump_pngs(hdict, label):
        for name, hh in hdict.items():
            outpng = os.path.join(outdir, "png", f"{hh.name}.png")
            xlabel = name.replace("_", " ")
            hh.save_png(outpng, xlabel=xlabel)

    dump_pngs(h_all, "ALL")
    dump_pngs(h_f5, "FID5")
    dump_pngs(h_f4, "FID4")

    # Write summary
    summary_path = os.path.join(outdir, f"{tag}_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"File: {lhe_path}\n")
        f.write(f"Tag: {tag}\n")
        f.write(f"Beams (pdg1,pdg2,E1,E2,NPRUP): {beams}\n")
        f.write(f"MGGenerationInfo: nevt={nevt_mg}, xsec_total(pb)={xsec_mg}\n")
        f.write(f"Init sum xsec(pb)={xsec_init}\n")
        f.write(f"Using total xsec(pb)={xsec_total}\n\n")
        f.write(f"Events read: {nev_read}\n")
        f.write(f"Events reconstructed: {nev_reco}\n\n")
        f.write(f"Acceptance A(|eta|<5) = {A5:.12g}\n")
        f.write(f"Sigma_fid(|eta|<5) per-proc = {sigma_fid5:.12g} pb\n")
        f.write(f"Sigma_fid(|eta|<5) total*A  = {sigma_fid5_alt:.12g} pb\n\n")
        f.write(f"Acceptance A(|eta|<4) = {A4:.12g}\n")
        f.write(f"Sigma_fid(|eta|<4) per-proc = {sigma_fid4:.12g} pb\n")
        f.write(f"Sigma_fid(|eta|<4) total*A  = {sigma_fid4_alt:.12g} pb\n\n")
        f.write("Per-subprocess XSECUP:\n")
        for i, x in xsecs_proc.items():
            f.write(f"  proc {i}: {x:.12g} pb\n")

    print(f"\n=== {tag} ===")
    print(f"Total xsec = {xsec_total} pb")
    print(f"A(|eta|<5) = {A5:.6g}  -> sigma_fid5 = {sigma_fid5:.12g} pb")
    print(f"A(|eta|<4) = {A4:.6g}  -> sigma_fid4 = {sigma_fid4:.12g} pb")
    print(f"Outputs: {summary_path} | {npz_path} | {os.path.join(outdir,'png')}")

    return {
        "tag": tag,
        "xsec_total_pb": xsec_total,
        "A5": A5, "sigma_fid5": sigma_fid5,
        "A4": A4, "sigma_fid4": sigma_fid4,
        "summary": summary_path,
        "npz": npz_path,
        "png_dir": os.path.join(outdir, "png"),
    }


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--lhec", default=LHEC_LHE)
    ap.add_argument("--lhmuc", default=LHMUC_LHE)
    ap.add_argument("--outdir", default="ttbar_genlevel_outputs_noROOT")
    ap.add_argument("--maxev", type=int, default=None, help="quick test: only process N events")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    analyze_sample(args.lhec, os.path.join(args.outdir, "LHeC"), "LHeC", max_events=args.maxev)
    analyze_sample(args.lhmuc, os.path.join(args.outdir, "LHmuC"), "LHmuC", max_events=args.maxev)


if __name__ == "__main__":
    main()
