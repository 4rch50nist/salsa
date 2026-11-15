import numpy as np

G, WIN = 40, 250
WORDS = "zero one two three four five six seven eight nine ten eleven twelve".split()

def curve(stall, pos=None):
    w = np.array([1.0 / float(l.split()[1]) for l in open(stall) if l.startswith("W ")])
    rows = []
    for g in range(len(w) // G):
        c = w[g * G:(g + 1) * G]
        rows.append([(g + 0.5) * G * WIN / 1e6, np.median(c), np.percentile(c, 10), np.percentile(c, 90)])
    d = np.array(rows)
    if pos is not None:
        p = np.loadtxt(pos, dtype=np.int64, ndmin=1)
        d[:, 0] = p[np.minimum(2 * np.arange(len(d)) + 1, len(p) - 1)] / 1e6
    return d

def table(d):
    return "x med lo hi\n" + "\n".join(f"{x:.4f} {m:.4f} {l:.4f} {h:.4f}" for x, m, l, h in d) + "\n}"

def series(d, tag, band, line):
    return (f"  \\addplot[draw=none, name path={tag}lo, forget plot] table[x=x, y=lo]{{\n{table(d)};\n"
            f"  \\addplot[draw=none, name path={tag}hi, forget plot] table[x=x, y=hi]{{\n{table(d)};\n"
            f"  \\addplot[{band}, fill opacity=0.35, draw=none, forget plot] fill between[of={tag}lo and {tag}hi];\n"
            f"  \\addplot[{line}, forget plot] table[x=x, y=med]{{\n{table(d)};\n")

ns = int(open("genome_meta.txt").read().split()[0])
nw = WORDS[ns] if ns < len(WORDS) else str(ns)
A = curve("genome_stall.txt")
F = curve("genome_focus_stall.txt", "genome_focus_pos.txt")
xmax = np.ceil(max(A[-1, 0], F[-1, 0]) * 20) / 20
ymax = max(1.65, np.ceil((max(A[:, 3].max(), F[:, 3].max()) + 0.05) * 20) / 20)
ticks = ",".join(f"{t:.1f}" for t in np.arange(1.0, ymax + 1e-9, 0.2))

tex = (r"""\begin{figure}[!t]
\centering
\begin{tikzpicture}[trim axis right]
\begin{axis}[width=\columnwidth, height=0.5\columnwidth, axis lines=box, axis line style={black, line width=0.5pt}, tick align=inside, tick style={black!70},
  grid=both, minor x tick num=1, minor y tick num=1, scaled x ticks=false, xticklabel style={/pgf/number format/fixed},
  major grid style={black!40, dotted, line width=0.5pt}, minor grid style={black!18, dotted, line width=0.4pt},
  tick label style={font=\scriptsize}, label style={font=\footnotesize}, title style={font=\footnotesize},
  ylabel={Cycles per Operation}, ylabel style={yshift=-2pt},
  legend style={font=\scriptsize, draw=black!60, fill=white, fill opacity=0.9, text opacity=1, legend cell align=left},
""" + f"  xmin=0, xmax={xmax:g}, ymin=0.95, ymax={ymax:g}, ytick={{{ticks}}},\n"
    + r"""  xlabel={Operations (Millions)}, title={Genome coverage},
  legend style={at={(0.98,0.97)}, anchor=north east, legend columns=2, column sep=5pt}]
""" + f"  \\draw[black!60, line width=0.6pt, loosely dashed] (axis cs:0,1) -- (axis cs:{xmax:g},1);\n"
    + series(A, "ga", "fill=black!18", "black, line width=1.0pt")
    + series(F, "gf", "fill=black!10", "black!60, line width=1.0pt, densely dashed")
    + f"  \\addlegendimage{{black, line width=1.0pt}}\\addlegendentry{{All {nw} samples}}\n"
    + "  \\addlegendimage{black!60, line width=1.0pt, densely dashed}\\addlegendentry{Busiest sample alone}\n"
    + r"""\end{axis}
\end{tikzpicture}
""" + f"\\caption{{Cycles per operation of a pipeline that stalls whenever an operation needs a write still in flight, over SARS-CoV-2 amplicon coverage for all {nw} samples of a multiplexed ARTIC nanopore run and for the busiest sample alone. The dashed line marks one cycle per operation; lines are medians, and bands span the 10th to 90th percentiles.}}\n"
    + "\\label{fig:genome_hazards}\n\\end{figure}\n")
open("genome_fig.tex", "w").write(tex)
print("wrote genome_fig.tex")
