"""Wong colorblind-friendly palette matching Makie defaults."""

# Wong (2011) palette — same cycle Makie uses by default.
WONG = [
    (0.0, 0.4470, 0.7410),      # blue
    (0.8500, 0.3250, 0.0980),   # orange
    (0.9290, 0.6940, 0.1250),   # yellow
    (0.4940, 0.1840, 0.5560),   # purple
    (0.4660, 0.6740, 0.1880),   # green
    (0.3010, 0.7450, 0.9330),   # sky blue
    (0.6350, 0.0780, 0.1840),   # dark red
]


def particle_color(i: int) -> tuple[float, float, float]:
    """Return RGB tuple for the i-th particle (1-indexed, Makie parity)."""
    return WONG[(i - 1) % len(WONG)]


def particle_color_rgba(i: int, alpha: float = 1.0) -> tuple[float, float, float, float]:
    r, g, b = particle_color(i)
    return (r, g, b, alpha)


def particle_color_255(i: int) -> tuple[int, int, int]:
    r, g, b = particle_color(i)
    return (int(r * 255), int(g * 255), int(b * 255))
