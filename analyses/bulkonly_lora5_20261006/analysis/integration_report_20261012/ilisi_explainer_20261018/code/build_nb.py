import nbformat as nbf
C=[]
C.append(nbf.v4.new_markdown_cell("""# iLISI 검정 설명 그림 (2026-10-18)
RF 설명 슬라이드와 같은 흐름: ① 개념 (toy 예시, 식, 읽는 법) → 파이프라인 → 기준값 척도, ② 대조군 설계와 실제 결과 (sepsis 공동 임베딩 UMAP, pre→post Harmony, 환자 1:1).
입력 `inputs/` = 통합 보고서 `results/check/` (calibration_summary.csv, hvgpca_sensitivity_summary.csv, *_umap_*.npz)와 `patient_pair_ilisi_20261018/results/`.
LISI 계산은 통합 보고서 `code/metrics.py`의 `_lisi_from_knn`와 같은 식 (perplexity 30, k = 90)."""))
C.append(nbf.v4.new_code_cell(r'''import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, os
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from sklearn.neighbors import NearestNeighbors
BASE = os.environ.get("IL_BASE", "."); IN = f"{BASE}/inputs"; FIG = f"{BASE}/figs"; RES = f"{BASE}/results"
for d in (FIG, RES): os.makedirs(d, exist_ok=True)
mpl.rcParams.update({"font.family": ["Noto Sans CJK KR", "DejaVu Sans"], "font.size": 10, "axes.unicode_minus": False,
    "axes.edgecolor": "#B8BDC7", "axes.labelcolor": "#4A5162", "xtick.color": "#4A5162", "ytick.color": "#4A5162",
    "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#FFFFFF", "savefig.facecolor": "#FFFFFF"})
INK, INK2, MUTED, GRID = "#1F2430", "#4A5162", "#8A90A0", "#E4E6EB"
REF = "#B8BDC7"
CMP = {"realU": ("L1", "같은 연구 다른 표본 real", "#1baf7a"), "realT1": ("L2", "같은 표본, depth만 맞춤", "#008300"),
       "other": ("L3", "다른 연구 real", "#6250d6"), "realT2": ("L4", "PBISC 형식으로 바꾼 real", "#e87ba4"),
       "k1000": ("L5", "PBISC K = 1,000", "#2a78d6"), "k4000": ("L6", "PBISC K = 4,000", "#eb6834")}

def lisi_from_knn(d, idx, labels, perplexity=30, tol=1e-5, iters=50):
    """identical to integration_report_20261012/code/metrics.py::_lisi_from_knn"""
    lab = pd.Categorical(labels).codes; L = lab.max() + 1; n = len(d); out = np.zeros(n); logU = np.log(perplexity)
    for i in range(n):
        di = d[i] ** 2; beta = 1.0; lo, hi = -np.inf, np.inf
        for _ in range(iters):
            P = np.exp(-di * beta); s = P.sum()
            if s == 0: H = 0.0; P = np.ones_like(di) / len(di)
            else: H = np.log(s) + beta * (di * P).sum() / s; P = P / s
            diff = H - logU
            if abs(diff) < tol: break
            if diff > 0: lo = beta; beta = beta * 2 if hi == np.inf else (beta + hi) / 2
            else: hi = beta; beta = beta / 2 if lo == -np.inf else (beta + lo) / 2
        w = np.bincount(lab[idx[i]], weights=P, minlength=L); out[i] = 1.0 / (w ** 2).sum()
    return out'''))
C.append(nbf.v4.new_markdown_cell("## 1. 실제 값 정리"))
C.append(nbf.v4.new_code_cell(r'''cal = pd.read_csv(f"{IN}/calibration_summary.csv")
KEY = {"Real, held-out samples (same study)": "realU", "Real, depth-matched": "realT1", "Real, PBISC-like read-out": "realT2",
       "Real, other study (COVID)": "other", "Real, other study (sepsis)": "other", "Generated K = 1,000": "k1000", "Generated K = 4,000": "k4000"}
cal["cmp"] = cal.source.map(KEY); val = cal.dropna(subset=["cmp"])[["cohort", "cmp", "iLISI_pre", "iLISI_post"]]
val.to_csv(f"{RES}/ilisi_pooled_values.csv", index=False)
pp = pd.read_csv(f"{IN}/pair_ilisi_per_repeat.csv"); pg = pd.read_csv(f"{IN}/pair_ilisi_real_vs_generated_per_repeat.csv")
pair = pd.concat([pp, pg], ignore_index=True); pair.to_csv(f"{RES}/ilisi_patient_pair_values.csv", index=False)
val.pivot(index="cmp", columns="cohort", values=["iLISI_pre", "iLISI_post"]).round(3)'''))
C.append(nbf.v4.new_markdown_cell("## 2. Toy 예시: 같은 LISI 식으로 계산"))
C.append(nbf.v4.new_code_cell(r'''rng = np.random.default_rng(7); n = 400
def cloud(shift):
    a = rng.normal(0, 1, (n, 2)); b = rng.normal(0, 1, (n, 2)) + shift
    return np.vstack([a, b]), np.r_[np.zeros(n, int), np.ones(n, int)]
TOY = {}
for name, sh in (("sep", [4.2, 0]), ("part", [1.3, 0.4]), ("mix", [0, 0])):
    X, s = cloud(np.array(sh)); nn = NearestNeighbors(n_neighbors=90).fit(X); d, idx = nn.kneighbors(X)
    li = lisi_from_knn(d, idx, s); TOY[name] = (X, s, d, idx, li)
    print(name, "median iLISI", round(float(np.median(li)), 2))'''))
C.append(nbf.v4.new_code_cell(r'''def toy_panel(ax, name, title):
    X, s, d, idx, li = TOY[name]
    # focal cell = the real cell whose iLISI is closest to the median
    cand = np.where(s == 0)[0]; f = cand[np.argmin(np.abs(li[cand] - np.median(li)))]
    nb = idx[f]; r = d[f].max()
    ax.scatter(*X[s == 0].T, s=7, color=REF, lw=0); ax.scatter(*X[s == 1].T, s=7, color="#2a78d6", alpha=0.55, lw=0)
    ax.scatter(*X[nb].T, s=13, facecolor="none", edgecolor=INK, lw=0.5)
    ax.add_patch(Circle(X[f], r, fill=False, ec=INK, lw=1.0, ls="--"))
    ax.scatter(*X[f], s=70, marker="*", color="#e34948", zorder=5, edgecolor="white", lw=0.6)
    pg = (s[nb] == 1).mean()
    ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_xlim(-3.6, 7.6); ax.set_ylim(-3.6, 3.9); ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=10.5, color=INK, fontweight="bold", pad=4)
    ax.text(0.5, -0.03, f"★ 세포의 이웃 90개 중 비교 세포 {100*pg:.0f}%", transform=ax.transAxes, ha="center", va="top", fontsize=9, color=INK2)
    ax.text(0.5, -0.17, f"중앙값 iLISI = {np.median(li):.2f}", transform=ax.transAxes, ha="center", va="top", fontsize=11.5, color="#2a78d6", fontweight="bold")

def box(ax, x, y, w, h, title, body, ec="#C9CED8", fc="#FFFFFF", tsize=10.5, bsize=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.015", fc=fc, ec=ec, lw=1.3, transform=ax.transAxes, clip_on=False))
    ax.text(x + 0.02, y + h - 0.05, title, fontsize=tsize, fontweight="bold", color=INK, va="top", transform=ax.transAxes)
    ax.text(x + 0.02, y + h - 0.13, body, fontsize=bsize, color=INK2, va="top", transform=ax.transAxes, linespacing=1.45)'''))
C.append(nbf.v4.new_markdown_cell("## 3. 그림 ① 개념 · 파이프라인 · 읽는 법"))
C.append(nbf.v4.new_code_cell(r'''fig = plt.figure(figsize=(16, 9), dpi=100)
fig.text(0.045, 0.955, "iLISI · 개념: 각 세포의 이웃에 real과 생성 세포가 얼마나 섞였나", fontsize=19, fontweight="bold", color=INK)
fig.text(0.045, 0.918, "LISI (Korsunsky et al. Nat Methods 2019, Harmony 논문) · 두 출처일 때 1 = 이웃이 한쪽뿐, 2 = 반반 · 생성 모델 논문에서는 GRouNdGAN (Zinati 2024)이 사용",
         fontsize=10.5, color=INK2)
# A toy
gsA = fig.add_gridspec(1, 3, left=0.045, right=0.60, top=0.83, bottom=0.57, wspace=0.10)
fig.text(0.045, 0.877, "A  같은 식을 toy 데이터에 적용 (출처당 400 세포)", fontsize=12.5, fontweight="bold", color=INK)
for j, (nm, tt) in enumerate((("sep", "완전 분리"), ("part", "일부 겹침"), ("mix", "같은 분포"))):
    toy_panel(fig.add_subplot(gsA[0, j]), nm, tt)
fig.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", color=REF, label="기준 real"),
                    mpl.lines.Line2D([], [], marker="o", ls="", color="#2a78d6", label="비교 세포"),
                    mpl.lines.Line2D([], [], marker="*", ls="", color="#e34948", ms=11, label="한 세포"),
                    mpl.lines.Line2D([], [], ls="--", color=INK, label="그 세포의 이웃 90개")],
           loc="upper left", bbox_to_anchor=(0.045, 0.862), ncol=4, frameon=False, fontsize=9)
# B formula + reading
axB = fig.add_axes([0.63, 0.47, 0.335, 0.40]); axB.axis("off")
box(axB, 0.0, 0.62, 1.0, 0.38, "식 (세포 한 개마다)",
    "이웃 k = 90개를 거리에 따라 가중 (perplexity 30)\n출처별 가중 비율 p_real, p_gen\niLISI = 1 / (p_real² + p_gen²)   → 반반이면 1/(0.25+0.25) = 2", ec="#2a78d6", fc="#EEF1F8")
box(axB, 0.0, 0.31, 1.0, 0.28, "값 읽는 법", "세포별 iLISI의 중앙값을 보고 (miLISI)\n1.0 = 이웃이 모두 같은 출처 · 2.0 = 완전히 섞임")
box(axB, 0.0, 0.0, 1.0, 0.28, "왜 real–real 대조가 필요한가",
    "환자·연구·기기·세포 수에 따라 상한이 2보다 낮음\n같은 파이프라인으로 real끼리 잰 값이 실제 기준선")
# C pipeline
axC = fig.add_axes([0.045, 0.235, 0.92, 0.18]); axC.axis("off")
fig.text(0.045, 0.43, "B  우리 파이프라인 (통합 보고서 c_joint.py · e_checks.py)", fontsize=12.5, fontweight="bold", color=INK)
steps = [("기준 real", "표본 절반 A\n표본당 400 (sepsis)\n/ 300 (COVID) 세포"), ("비교 세포 X", "real B, 다른 연구,\nPBISC 형식 real,\nPBISC K1k / K4k"),
         ("정규화", "공유 패널 5,675 유전자\nlog1p CP10k\n(패널 총합 기준)"), ("HVG · PCA", "HVG 2,000 (seurat_v3)\nPCA 30을 기준에 fit,\nX는 투영"),
         ("Harmony", "batch = 출처\n(기준 vs X, 쌍마다)\n전 = pre, 후 = post"), ("iLISI", "k = 90, perplexity 30\n세포 12,000개 표본\n중앙값 보고")]
w = 0.148; gap = (1 - 6 * w) / 5
for i, (t, b) in enumerate(steps):
    x = 0.004 + i * (w + gap) * 0.99; hl = t in ("Harmony", "iLISI")
    axC.add_patch(FancyBboxPatch((x, 0.05), w, 0.9, boxstyle="round,pad=0.006,rounding_size=0.03", fc="#EEF1F8" if hl else "#FFFFFF",
                                 ec="#2a78d6" if hl else "#C9CED8", lw=1.3, transform=axC.transAxes, clip_on=False))
    axC.text(x + w / 2, 0.82, t, ha="center", va="top", fontsize=11, fontweight="bold", color=INK, transform=axC.transAxes)
    axC.text(x + w / 2, 0.58, b, ha="center", va="top", fontsize=8.8, color=INK2, transform=axC.transAxes, linespacing=1.4)
    if i < 5: axC.add_patch(FancyArrowPatch((x + w + 0.004, 0.5), (x + w + gap - 0.004, 0.5), arrowstyle="-|>", mutation_scale=13, color="#6A7080", transform=axC.transAxes))
# D scale
axD = fig.add_axes([0.065, 0.06, 0.87, 0.10])
fig.text(0.045, 0.19, "C  기준값 척도 (Harmony 후 중앙값 · sepsis 통합 보고서 V0, GRouNdGAN은 논문 값)", fontsize=12.5, fontweight="bold", color=INK)
axD.set_xlim(0.95, 2.05); axD.set_ylim(0, 1); axD.axis("off")
axD.add_patch(mpl.patches.FancyBboxPatch((1.0, 0.42), 1.0, 0.16, boxstyle="round,pad=0,rounding_size=0.02", fc="#EEF1F5", ec="none"))
grad = np.linspace(0, 1, 256)[None, :]; axD.imshow(grad, extent=(1.0, 2.0, 0.42, 0.58), aspect="auto", cmap=mpl.colors.LinearSegmentedColormap.from_list("g", ["#DCE6F5", "#2a78d6"]), alpha=0.55)
axD.text(1.0, 0.2, "1.0 · 이웃이 한 출처", ha="left", fontsize=9.5, color=INK2); axD.text(2.0, 0.2, "2.0 · 완전히 섞임", ha="right", fontsize=9.5, color=INK2)
s = val[val.cohort == "sepsis"].set_index("cmp").iLISI_post
marks = [(1.94, "GRouNdGAN 생성 vs real\n(같은 데이터로 학습·평가)", "#8A90A0", 1), (s["realU"], "같은 연구 다른 표본 real", CMP["realU"][2], 0),
         (s["other"], "다른 연구 real", CMP["other"][2], 1), (s["realT2"], "PBISC 형식 real", CMP["realT2"][2], 0),
         ((s["k1000"] + s["k4000"]) / 2, f"PBISC {s['k1000']:.2f}–{s['k4000']:.2f}", CMP["k1000"][2], 1)]
for v, lab, col, up in marks:
    axD.plot([v, v], [0.40, 0.60], color=col, lw=3)
    axD.text(v, 0.68 if up else 0.03, f"{v:.2f}  {lab}" if "PBISC " not in lab or "형식" in lab else lab, ha="center", va="bottom" if up else "bottom",
             fontsize=9, color=INK, fontweight="bold" if "PBISC" in lab and "형식" not in lab else "normal")
fig.savefig(f"{FIG}/ilisi_1_concept.png", dpi=300); plt.show()'''))
C.append(nbf.v4.new_markdown_cell("## 4. 그림 ② 설계와 결과"))
C.append(nbf.v4.new_code_cell(r'''def umap_panel(ax, coh, key, title, v):
    z = np.load(f"{IN}/{coh}_umap_{key}.npz", allow_pickle=True); U, s = z["umap"], z["src"]
    o = np.random.default_rng(0).permutation(len(U)); U, s = U[o], s[o]
    col = CMP["other" if key in ("covid", "sepsis") else key][2]
    ax.scatter(*U.T, s=0.7, c=np.where(s == 0, REF, col), lw=0, rasterized=True)   # shuffled order: neither source is drawn on top
    ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(False) for sp in ax.spines.values()]
    ax.set_title(title, fontsize=10.5, color=INK, fontweight="bold", pad=3)
    ax.text(0.5, -0.02, f"iLISI {v:.2f}", transform=ax.transAxes, ha="center", va="top", fontsize=12, color=col, fontweight="bold")

def dot_panel(ax, coh, show_y=True):
    order = ["realU", "realT1", "other", "realT2", "k1000", "k4000"]; y = np.arange(len(order))[::-1]
    d = val[val.cohort == coh].set_index("cmp")
    for yy, k in zip(y, order):
        a, b = d.loc[k, "iLISI_pre"], d.loc[k, "iLISI_post"]; col = CMP[k][2]
        ax.add_patch(FancyArrowPatch((a, yy), (b, yy), arrowstyle="-|>", mutation_scale=10, color=col, lw=1.6, shrinkA=4, shrinkB=5))
        ax.scatter(a, yy, s=48, facecolor="white", edgecolor=col, lw=1.6, zorder=3); ax.scatter(b, yy, s=60, color=col, zorder=4)
        ax.text(b + 0.025, yy + 0.18, f"{b:.2f}", fontsize=9, color=INK, va="center")
    ax.set_yticks(y, [f"{CMP[k][0]}  {CMP[k][1]}" for k in order] if show_y else [], fontsize=9.5)
    ax.set_xlim(0.95, 1.95); ax.set_xticks([1.0, 1.25, 1.5, 1.75]); ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.set_ylim(-0.6, len(order) - 0.4); ax.set_xlabel("iLISI 중앙값 (기준 real vs X)")
    ax.set_title({"sepsis": "Sepsis (Ye 2026)", "covid": "COVID (Stephenson 2025, 187 표본)"}[coh], loc="left", fontsize=11, fontweight="bold", color=INK)

fig = plt.figure(figsize=(16, 9), dpi=100)
fig.text(0.045, 0.955, "iLISI · 설계와 결과: real끼리의 값과 나란히 놓고 읽는다", fontsize=19, fontweight="bold", color=INK)
fig.text(0.045, 0.918, "기준 real = 표본 절반 A · 비교 X는 L1–L6 · 쌍마다 따로 Harmony (batch = 출처) · 회색 = 기준 real", fontsize=10.5, color=INK2)
fig.text(0.045, 0.875, "A  Harmony 후 공동 UMAP (sepsis, 출처당 15,000 세포)", fontsize=12.5, fontweight="bold", color=INK)
gsU = fig.add_gridspec(1, 4, left=0.045, right=0.975, top=0.85, bottom=0.555, wspace=0.06)
s = val[val.cohort == "sepsis"].set_index("cmp").iLISI_post
for j, (key, tt, v) in enumerate((("realU", "L1  같은 연구 다른 표본 real", s["realU"]), ("covid", "L3  다른 연구 real (COVID)", s["other"]),
                                  ("realT2", "L4  PBISC 형식 real", s["realT2"]), ("k1000", "L5  PBISC K = 1,000", s["k1000"]))):
    umap_panel(fig.add_subplot(gsU[0, j]), "sepsis", key, tt, v)
fig.text(0.045, 0.475, "B  pooled 비교: Harmony 전 ○ → 후 ●", fontsize=12.5, fontweight="bold", color=INK)
gsD = fig.add_gridspec(1, 2, left=0.19, right=0.62, top=0.43, bottom=0.12, wspace=0.08)
dot_panel(fig.add_subplot(gsD[0, 0]), "sepsis", True); dot_panel(fig.add_subplot(gsD[0, 1]), "covid", False)
# C patient 1:1
fig.text(0.665, 0.475, "C  환자 1:1 (sepsis, 400 + 400 세포, 5회 반복)", fontsize=12.5, fontweight="bold", color=INK)
axP = fig.add_axes([0.80, 0.12, 0.175, 0.29])
rows = [("same patient split", "같은 환자 반반 (상한)", "#4A4F5C"), ("all patients", "real 환자 vs real 환자", "#1baf7a"),
        ("real patient vs k1000 sample", "real 환자 vs PBISC K1k 표본", CMP["k1000"][2]), ("real patient vs k4000 sample", "real 환자 vs PBISC K4k 표본", CMP["k4000"][2])]
y = np.arange(len(rows))[::-1]
for yy, (k, lab, col) in zip(y, rows):
    d = pair[pair.set == k]; a, b = d.iLISI_pre.mean(), d.iLISI_post.mean()
    axP.add_patch(FancyArrowPatch((a, yy), (b, yy), arrowstyle="-|>", mutation_scale=10, color=col, lw=1.6, shrinkA=4, shrinkB=5))
    axP.scatter(a, yy, s=48, facecolor="white", edgecolor=col, lw=1.6, zorder=3)
    axP.scatter(d.iLISI_post, np.full(len(d), yy) + np.linspace(-0.12, 0.12, len(d)), s=12, color=col, alpha=0.6, zorder=4)
    axP.scatter(b, yy, s=60, color=col, zorder=5); axP.text(b + 0.02, yy + 0.25, f"{b:.2f}", fontsize=9, color=INK)
axP.set_yticks(y, [r[1] for r in rows], fontsize=9.5); axP.tick_params(axis="y", length=0); axP.spines["left"].set_visible(False); axP.set_xlim(1.1, 1.9); axP.set_xticks([1.2, 1.4, 1.6, 1.8])
axP.set_ylim(-0.6, len(rows) - 0.4); axP.xaxis.grid(True, color=GRID); axP.set_axisbelow(True); axP.set_xlabel("iLISI (쌍별 중앙값의 평균)")
fig.text(0.665, 0.452, "400 + 400 세포 척도 → B (pooled)와 직접 비교하지 않음 · ○ 전 → ● 후, 작은 점 = 반복", fontsize=9, color=MUTED)
fig.text(0.045, 0.035, "읽기: Harmony 후 PBISC (sepsis 1.14–1.15, COVID 1.06)는 다른 연구 real (1.51)보다 낮고, 형식만 바꾼 real (1.44)과의 차이가 모델 잔차. "
         "환자 1:1 척도에서는 PBISC 표본이 real 환자쌍과 같은 수준 (1.47–1.48 vs 1.48).", fontsize=10, color=INK)
fig.savefig(f"{FIG}/ilisi_2_results.png", dpi=300); plt.show()'''))
nb=nbf.v4.new_notebook(); nb.cells=C; nb.metadata["kernelspec"]={"name":"python3","display_name":"Python 3","language":"python"}
nbf.write(nb,'/home/claude/ilisi/ilisi_explainer.ipynb'); print('ok')
