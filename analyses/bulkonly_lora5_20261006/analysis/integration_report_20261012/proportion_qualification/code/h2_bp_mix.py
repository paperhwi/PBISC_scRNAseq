"""Rebuild the BayesPrism+cal25 mixture with a per-sample calibration (25 LoRA training donors minus the sample itself),
matching the WNNLS+cal25 rule; the reference is unchanged."""
exec(open("h_equal_info.py").read().split("cats, sig, mk = markers_sig(ref)")[0])
R = pd.read_csv(f"{O}/bp/eq25cal_ref.csv.gz", index_col=0, nrows=1); genes = list(R.columns)
pos = {s_: i for i, s_ in enumerate(sym[gi])}; g = np.array([pos[x] for x in genes])
rows = {}
for s in ids:
    f, n = calib([x for x in D25 if x != s], g); rows[s] = np.rint(Bk[s].to_numpy()[g] * f)
mix = pd.DataFrame(rows, index=genes).T.astype(int); mix.to_csv(f"{O}/bp/eq25cal_mix.csv.gz")
m = json.load(open(f"{E}/equal_info_meta.json")); m["bayesprism_cal_donors"] = "25 LoRA training donors minus the sample itself (21-22 usable)"; json.dump(m, open(f"{E}/equal_info_meta.json", "w"), indent=1)
print(mix.shape)
