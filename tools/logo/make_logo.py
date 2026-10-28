"""Generate the kdensity hex logo (text still as <text>; render.sh outlines it).

The curves are exact sections u -> c(u | v; rho) of the Gaussian copula
density, the kernel of Jones & Henderson's gcopula estimator, with
rho = 1 - h^2. The colours come from the 2018 kDensity wordmark:
blue #00a2e8 and red #ed1c24, with the red "d" as a nod to its red "D".
"""

import math
import os
import sys
from statistics import NormalDist

N = NormalDist()

BLUE = "#00a2e8"
RED = "#ed1c24"
NAVY = "#0b3f63"
WHITE = "#ffffff"

R = 300.0                      # circumradius of the pointy-top hexagon
W = R * math.sqrt(3)           # width, about 519.6
H = 2 * R                      # height, 600
CX, CY = W / 2, R
BORDER = 20


def hexagon(r):
    pts = []
    for k in range(6):
        t = math.radians(90 + 60 * k)
        pts.append((CX + r * math.cos(t), CY - r * math.sin(t)))
    return "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"


def inner_r(border):
    return R - border / math.cos(math.radians(30))


def section(v, h, n=600):
    """Samples (u, c(u | v)) of a Gaussian-copula section on [0, 1]."""
    rho = 1 - h * h
    b = N.inv_cdf(v)
    den = 2 * (1 - rho * rho)
    norm = 1 / math.sqrt(1 - rho * rho)
    out = [(0.0, 0.0)]
    for i in range(1, n):
        # Chebyshev spacing: dense near the edges, where the sections bend fastest.
        u = 0.5 - 0.5 * math.cos(math.pi * i / n)
        a = N.inv_cdf(u)
        out.append((u, norm * math.exp(-(rho * rho * (a * a + b * b) - 2 * rho * a * b) / den)))
    out.append((1.0, 0.0))
    return out


def to_path(samples, x0, x1, base, scale, close=False):
    pts = [(x0 + u * (x1 - x0), base - y * scale) for u, y in samples]
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    if close:
        d += " Z"
    return d


def wordmark(y, size):
    font = (f"font-family:'Comfortaa';font-weight:bold;font-size:{size}px;"
            f"letter-spacing:1px;text-anchor:middle;stroke-width:3px;stroke-linejoin:round")

    def paint(col):
        return f"fill:{col};stroke:{col}"

    return (f'<text x="{CX:.2f}" y="{y}" style="{font};{paint(BLUE)}">'
            f'k<tspan style="{paint(RED)}">d</tspan>ensity</text>')


def logo():
    h = 0.45
    vs = [0.12, 0.27, 0.5, 0.73, 0.88]
    hi = 1                                 # the section drawn in red
    x0, x1, base = 62.0, W - 62.0, 352.0
    secs = [section(v, h) for v in vs]
    peak = max(max(y for _, y in s) for s in secs)
    scale = 215.0 / peak

    lines = []
    for i, s in enumerate(secs):
        if i != hi:
            lines.append(f'  <path d="{to_path(s, x0, x1, base, scale)}" fill="none" '
                         f'stroke="{BLUE}" stroke-width="6" stroke-opacity="0.55" stroke-linejoin="round"/>')
    s = secs[hi]
    lines.append(f'  <path d="{to_path(s, x0, x1, base, scale, close=True)}" fill="{RED}" fill-opacity="0.12"/>')
    lines.append(f'  <path d="{to_path(s, x0, x1, base, scale)}" fill="none" stroke="{RED}" '
                 f'stroke-width="10" stroke-linejoin="round" stroke-linecap="round"/>')
    inner = hexagon(inner_r(BORDER))
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:.2f} {H:.0f}" width="{W:.0f}" height="{H:.0f}">
<title>kdensity</title>
<defs><clipPath id="in"><path d="{inner}"/></clipPath></defs>
<path d="{hexagon(R)}" fill="{NAVY}"/>
<path d="{inner}" fill="{WHITE}"/>
<g clip-path="url(#in)">
  <line x1="40" y1="{base}" x2="{W - 40:.1f}" y2="{base}" stroke="{NAVY}" stroke-width="4" stroke-linecap="round"/>
{chr(10).join(lines)}
</g>
{wordmark(440, 86)}
</svg>
"""


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "logo_text.svg")
    with open(out, "w") as fh:
        fh.write(logo())
