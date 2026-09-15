"""A classification as a page to browse, instead of a spreadsheet to scroll.

`classes.py --page` writes one of these.  It carries the whole classification
inline -- there is no server and nothing to fetch -- and answers the questions
`classview` answers, by filtering in the browser: by kind, by Coxeter
polynomial, by size, by quipu, by a relation string you remember.

Two things it draws that a table cannot, and they are the reason a page is worth
having at all:

* **the class' quipu**, from the `P^(m_0,...,m_r)_(k_0,...,k_{r+1})` parameters
  in its name -- the main string across, the cords hanging off their feet.  A
  quipu class *is* that tree, so the tree is the class' real name and the
  notation only spells it.
* **each LNA as its quiver**, `1 -> 2 -> ... -> n` with an arc over the span of
  every relation.  Whether two relations overlap, and by how much, is the thing
  the whole classification turns on, and it is invisible in `1;2;3|3;4;5;6` and
  obvious in a drawing.

`render` returns the page as the Artifact publish skeleton wants it -- no
doctype, head or body of its own.  `renderStandalone` wraps that for a file you
open in a browser.  Browsers imply the missing structure either way, so the
fragment opens locally too; the wrapper only adds the charset and viewport.
"""

import html
import json

from . import classview as cv

MEMBER_LIMIT = 6000
"""Rows of per-LNA detail to carry.  Catalan growth means n = 12 has 58786 of
them, which is a page nobody wants to load; past this the class summaries are
still complete and the drill-down says it was truncated."""


def _quipuParameters(name):
    """(k, m) out of 'P^(m_0,...,m_r)_(k_0,...,k_{r+1})', or None."""
    if not name.startswith("P^(") or ")_(" not in name or not name.endswith(")"):
        return None
    cords, _, gaps = name[3:-1].partition(")_(")
    try:
        m = [int(x) for x in cords.split(",")]
        k = [int(x) for x in gaps.split(",")]
    except ValueError:
        return None
    return (k, m) if len(k) == len(m) + 1 else None


def _relationSpans(relations):
    """[(start, end)] for '1;2;3|3;4;5;6', as vertex numbers."""
    spans = []
    for path in relations.split("|"):
        vertices = [v for v in path.split(";") if v]
        if len(vertices) >= 2:
            spans.append((int(vertices[0]), int(vertices[-1])))
    return spans


def collect(classification):
    """The classification as the page's data: classes, then members per class."""
    classes = classification.classes()
    carried = 0
    data = {
        "length": classification.length,
        "lnaCount": sum(c.size for c in classes),
        "classes": [],
        "kinds": {kind: list(counts) for kind, counts in classification.kindCounts().items()},
        "collisions": [
            {"polynomial": polynomial, "classes": [c.name for c in found]}
            for polynomial, found in sorted(classification.coxeterCollisions().items())
        ],
        "unsettled": len(classification.unsettled()),
        "memberLimit": MEMBER_LIMIT,
    }
    for mutationClass in classes:
        members = []
        if carried < MEMBER_LIMIT:
            for relations, path, _numbering in classification.members(mutationClass.name):
                members.append({
                    "relations": relations,
                    "spans": _relationSpans(relations),
                    "path": path,
                })
                carried += 1
                if carried >= MEMBER_LIMIT:
                    break
        data["classes"].append({
            "name": mutationClass.name,
            "size": mutationClass.size,
            "kind": mutationClass.kind,
            "form": mutationClass.hereditaryForm,
            "coxeter": mutationClass.coxeterPolynomial,
            "representative": mutationClass.representative,
            "repSpans": _relationSpans(mutationClass.representative),
            "quipu": _quipuParameters(mutationClass.hereditaryForm or mutationClass.name),
            "members": members,
            "membersComplete": len(members) == mutationClass.size,
        })
    return data


STYLE = """
<style>
  :root {
    --ground: #f4f6f6;
    --surface: #ffffff;
    --raised: #eef2f1;
    --ink: #16211f;
    --muted: #5f6e6c;
    --faint: #8a9997;
    --line: #dde3e2;
    --accent: #0f6c66;
    --accent-soft: #dbeeec;
    --warn: #8a5e0f;
    --warn-soft: #f6ecd8;
    --flag: #9b3b2f;
    --flag-soft: #f7e3df;
    --grid: #c9d3d1;
    --shadow: 0 1px 2px rgba(22, 33, 31, .06);
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      --ground: #0f1514;
      --surface: #171f1e;
      --raised: #1e2827;
      --ink: #e7edec;
      --muted: #97a6a4;
      --faint: #6f7d7b;
      --line: #263130;
      --accent: #63cac0;
      --accent-soft: #14312e;
      --warn: #d9ac58;
      --warn-soft: #302711;
      --flag: #e6907c;
      --flag-soft: #33201c;
      --grid: #38443f;
      --shadow: none;
    }
  }
  :root[data-theme="dark"] {
    --ground: #0f1514;
    --surface: #171f1e;
    --raised: #1e2827;
    --ink: #e7edec;
    --muted: #97a6a4;
    --faint: #6f7d7b;
    --line: #263130;
    --accent: #63cac0;
    --accent-soft: #14312e;
    --warn: #d9ac58;
    --warn-soft: #302711;
    --flag: #e6907c;
    --flag-soft: #33201c;
    --grid: #38443f;
    --shadow: none;
  }

  body {
    background: var(--ground);
    color: var(--ink);
    font-family: "IBM Plex Sans", system-ui, -apple-system, sans-serif;
    font-size: 15px;
    line-height: 1.5;
    -webkit-font-smoothing: antialiased;
  }
  .wrap { max-width: 68rem; margin: 0 auto; padding-inline: 20px; padding-block: 0 4rem; }

  h1 {
    font-family: Newsreader, Georgia, serif;
    font-weight: 500;
    font-size: clamp(1.9rem, 4.5vw, 2.7rem);
    line-height: 1.1;
    margin: 0;
    letter-spacing: -.01em;
    text-wrap: balance;
  }
  h1 .sub { color: var(--muted); font-style: italic; }
  h2 {
    font-family: Newsreader, Georgia, serif;
    font-weight: 500;
    font-size: 1.35rem;
    margin: 0;
    letter-spacing: -.005em;
  }
  .eyebrow {
    font-size: .7rem; font-weight: 600; letter-spacing: .12em;
    text-transform: uppercase; color: var(--faint);
  }
  code, .mono, .poly {
    font-family: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
    font-variant-ligatures: none;
  }

  header.page { padding-block: 2.5rem 1.25rem; display: flex; flex-direction: column; gap: .6rem; }
  .lede { color: var(--muted); max-width: 46rem; }

  /* --- the counts, as a rule of proportions rather than tiles --- */
  .tally { display: flex; flex-wrap: wrap; gap: .4rem 1.5rem; align-items: baseline; }
  .tally .n {
    font-family: "IBM Plex Mono", monospace; font-size: 1.05rem; font-weight: 600;
    font-variant-numeric: tabular-nums;
  }
  .tally .of { color: var(--muted); font-size: .85rem; }

  /* --- controls --- */
  .controls {
    position: sticky; top: env(safe-area-inset-top, 0px); z-index: 5;
    background: var(--ground);
    border-block: 1px solid var(--line);
    padding-block: .7rem;
    display: flex; flex-wrap: wrap; gap: .5rem .75rem; align-items: center;
  }
  .controls input[type="search"] {
    flex: 1 1 14rem; min-width: 0;
    font: inherit; color: inherit;
    background: var(--surface);
    border: 1px solid var(--line); border-radius: 3px;
    padding: .38rem .6rem;
  }
  .controls input[type="search"]::placeholder { color: var(--faint); }
  .chips { display: flex; flex-wrap: wrap; gap: .35rem; }
  .chip {
    font: inherit; font-size: .8rem; color: var(--muted);
    background: transparent; border: 1px solid var(--line); border-radius: 100px;
    padding: .18rem .62rem; cursor: pointer;
  }
  .chip[aria-pressed="true"] {
    color: var(--accent); border-color: var(--accent); background: var(--accent-soft);
    font-weight: 500;
  }
  .chip .count { color: var(--faint); font-variant-numeric: tabular-nums; }
  .chip[aria-pressed="true"] .count { color: var(--accent); }
  :focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }

  /* --- collisions --- */
  .collisions { margin-top: 2rem; border-left: 3px solid var(--warn); padding-left: 1rem; }
  .collisions p { margin: .3rem 0 .8rem; color: var(--muted); max-width: 46rem; }
  .collision { margin-bottom: .9rem; }
  .collision .poly { font-size: .78rem; color: var(--muted); display: block; overflow-x: auto; }
  .collision .pair { display: flex; flex-wrap: wrap; gap: .4rem; margin-top: .25rem; }
  .collision .pair code { font-size: .85rem; }

  /* --- the class list: bands, not cards --- */
  .list { margin-top: 1.75rem; }
  .band { border-bottom: 1px solid var(--line); }
  .band > summary {
    display: grid;
    grid-template-columns: 5.5rem 1fr auto;
    gap: .35rem 1rem;
    align-items: center;
    padding: .7rem .25rem;
    cursor: pointer;
    list-style: none;
  }
  .band > summary::-webkit-details-marker { display: none; }
  .band > summary:hover { background: var(--raised); }
  .band[open] > summary { background: var(--raised); }
  .tree { width: 5.5rem; height: 2.3rem; display: block; overflow: visible; }
  .ident { min-width: 0; }
  .ident .name { font-size: .95rem; font-weight: 500; word-break: break-all; }
  .ident .why { font-size: .78rem; color: var(--muted); word-break: break-word; }
  .ident .why code { font-size: .76rem; color: var(--ink); }
  .measure { display: flex; align-items: center; gap: .55rem; justify-self: end; }
  .bar { width: clamp(3rem, 14vw, 7rem); height: 6px; background: var(--line); border-radius: 100px; overflow: hidden; }
  .bar span { display: block; height: 100%; background: var(--accent); }
  .count-n {
    font-family: "IBM Plex Mono", monospace; font-variant-numeric: tabular-nums;
    font-size: .9rem; min-width: 3.2rem; text-align: right;
  }
  .pill {
    font-size: .68rem; letter-spacing: .04em; text-transform: uppercase;
    padding: .1rem .45rem; border-radius: 3px; white-space: nowrap;
  }
  .pill.quipu { color: var(--accent); background: var(--accent-soft); }
  .pill.canonical { color: var(--warn); background: var(--warn-soft); }
  .pill.flagged { color: var(--flag); background: var(--flag-soft); }
  .pill.other { color: var(--muted); background: var(--raised); }

  /* --- the drill-down --- */
  .detail { padding: .3rem .25rem 1.4rem; display: flex; flex-direction: column; gap: 1rem; }
  .facts { display: grid; grid-template-columns: max-content 1fr; gap: .2rem .9rem; font-size: .85rem; }
  .facts dt { color: var(--faint); }
  .facts dd { margin: 0; overflow-x: auto; }
  .members { display: flex; flex-direction: column; gap: .1rem; }
  .member {
    display: grid; grid-template-columns: 1fr auto; gap: .2rem .9rem;
    align-items: center; padding: .28rem 0; border-top: 1px dotted var(--line);
  }
  .member .quiver { display: block; overflow: visible; }
  .member .rel { font-size: .78rem; color: var(--muted); word-break: break-all; }
  .member .step {
    font-family: "IBM Plex Mono", monospace; font-size: .76rem; color: var(--faint);
    justify-self: end; white-space: nowrap;
  }
  .truncated { font-size: .8rem; color: var(--warn); }
  .empty { color: var(--muted); padding: 2rem .25rem; }

  footer.page {
    margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--line);
    font-size: .8rem; color: var(--faint);
  }
  footer.page code { font-size: .78rem; }

  @media (max-width: 34rem) {
    .band > summary { grid-template-columns: 4rem 1fr; }
    .measure { grid-column: 1 / -1; justify-self: start; }
    .tree { width: 4rem; }
  }
  @media (prefers-reduced-motion: reduce) { * { transition: none !important; } }
</style>
"""

SCRIPT = r"""
<script>
(function () {
  var data = window.__classification;
  var list = document.getElementById("list");
  var search = document.getElementById("search");
  var chips = Array.prototype.slice.call(document.querySelectorAll(".chip"));
  var kindFilter = null;

  var PILL = {
    "quipu": "quipu", "canonical type": "canonical",
    "not piecewise hereditary": "flagged", "contradictory": "flagged"
  };

  function svg(tag, attrs) {
    var node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (var key in attrs) { node.setAttribute(key, attrs[key]); }
    return node;
  }

  /* The quipu: main string across, cords hanging from their feet. */
  function drawQuipu(parameters) {
    var box = svg("svg", {"class": "tree", viewBox: "0 0 88 36",
                          role: "img", "aria-hidden": "true"});
    if (!parameters) {
      /* No quipu to draw, and a blank column reads as a drawing that failed
         rather than one there is nothing to draw.  Say so. */
      var mark = svg("text", {x: 44, y: 22, "text-anchor": "middle",
                              "font-size": 9, fill: "var(--faint)"});
      mark.textContent = "no tree";
      box.appendChild(mark);
      return box;
    }
    var k = parameters[0], m = parameters[1];
    var main = [], feet = [];
    for (var i = 0; i < k[0]; i++) { main.push(null); }
    for (var c = 0; c < m.length; c++) {
      feet.push({at: main.length, cord: m[c]});
      main.push("foot");
      for (var g = 0; g < k[c + 1]; g++) { main.push(null); }
    }
    /* One step has to fit the main string across and the deepest cord down,
       or a long cord draws straight out of the box. */
    var cordMax = Math.max.apply(null, m.concat([0]));
    var span = Math.max(main.length - 1, 1);
    var step = Math.min(9, 82 / span, cordMax ? 26 / cordMax : 9);
    var x0 = (88 - span * step) / 2;
    var y0 = (36 - cordMax * step) / 2;
    var line = function (x1, y1, x2, y2) {
      box.appendChild(svg("line", {x1: x1, y1: y1, x2: x2, y2: y2,
                                   stroke: "var(--grid)", "stroke-width": 1}));
    };
    var dot = function (x, y, r, fill) {
      box.appendChild(svg("circle", {cx: x, cy: y, r: r, fill: fill}));
    };
    for (var i = 0; i + 1 < main.length; i++) {
      line(x0 + i * step, y0, x0 + (i + 1) * step, y0);
    }
    for (var i = 0; i < main.length; i++) {
      var isFoot = main[i] === "foot";
      dot(x0 + i * step, y0, isFoot ? 2.1 : 1.5,
          isFoot ? "var(--accent)" : "var(--muted)");
    }
    feet.forEach(function (foot) {
      var x = x0 + foot.at * step;
      for (var j = 0; j < foot.cord; j++) {
        line(x, y0 + j * step, x, y0 + (j + 1) * step);
        dot(x, y0 + (j + 1) * step, 1.5, "var(--muted)");
      }
    });
    return box;
  }

  /* An LNA: 1 -> 2 -> ... -> n, with an arc over every relation's span. */
  function drawQuiver(length, spans) {
    var step = Math.min(13, 300 / Math.max(length - 1, 1));
    var width = (length - 1) * step + 8;
    var depth = spans.length ? 12 : 3;
    var box = svg("svg", {"class": "quiver", viewBox: "0 0 " + width + " " + (depth + 10),
                          width: width, height: depth + 10,
                          role: "img",
                          "aria-label": length + " vertices, " + spans.length + " relations"});
    var baseline = depth + 4, x0 = 4;
    for (var i = 0; i + 1 < length; i++) {
      box.appendChild(svg("line", {x1: x0 + i * step, y1: baseline,
                                   x2: x0 + (i + 1) * step, y2: baseline,
                                   stroke: "var(--grid)", "stroke-width": 1}));
    }
    for (var i = 0; i < length; i++) {
      box.appendChild(svg("circle", {cx: x0 + i * step, cy: baseline, r: 1.4,
                                     fill: "var(--muted)"}));
    }
    spans.forEach(function (pair) {
      var a = x0 + (pair[0] - 1) * step, b = x0 + (pair[1] - 1) * step;
      var lift = Math.min(depth, 3 + (b - a) * 0.34);
      box.appendChild(svg("path", {
        d: "M " + a + " " + baseline + " Q " + ((a + b) / 2) + " " + (baseline - lift * 2)
           + " " + b + " " + baseline,
        fill: "none", stroke: "var(--accent)", "stroke-width": 1.1
      }));
    });
    return box;
  }

  function renderMembers(mutationClass) {
    var holder = document.createElement("div");
    holder.className = "detail";

    var facts = document.createElement("dl");
    facts.className = "facts";
    var rows = [["Coxeter polynomial", mutationClass.coxeter]];
    if (mutationClass.form && mutationClass.form !== mutationClass.name) {
      rows.unshift(["Hereditary form", mutationClass.form]);
    }
    rows.push(["Representative", mutationClass.representative || "no relations"]);
    rows.forEach(function (pair) {
      var dt = document.createElement("dt");
      dt.textContent = pair[0];
      var dd = document.createElement("dd");
      var code = document.createElement("code");
      code.textContent = pair[1];
      dd.appendChild(code);
      facts.appendChild(dt);
      facts.appendChild(dd);
    });
    holder.appendChild(facts);

    var members = document.createElement("div");
    members.className = "members";
    mutationClass.members.forEach(function (member) {
      var row = document.createElement("div");
      row.className = "member";
      var left = document.createElement("div");
      left.appendChild(drawQuiver(data.length, member.spans));
      var rel = document.createElement("div");
      rel.className = "rel mono";
      rel.textContent = member.relations || "no relations";
      left.appendChild(rel);
      row.appendChild(left);
      var step = document.createElement("div");
      step.className = "step";
      step.textContent = member.path ? "mutate " + member.path.replace(/;/g, ", ")
                                     : "seeded";
      step.title = member.path ? "the mutations that reach this LNA from the class representative"
                               : "placed by the quipu theorem, with no search";
      row.appendChild(step);
      members.appendChild(row);
    });
    holder.appendChild(members);

    if (!mutationClass.membersComplete) {
      var note = document.createElement("p");
      note.className = "truncated";
      note.textContent = "Showing " + mutationClass.members.length + " of "
        + mutationClass.size + " -- the page carries at most " + data.memberLimit
        + " rows of detail. Use classes.py --members for the rest.";
      holder.appendChild(note);
    }
    return holder;
  }

  function band(mutationClass) {
    var details = document.createElement("details");
    details.className = "band";
    var summary = document.createElement("summary");

    summary.appendChild(drawQuipu(mutationClass.quipu));

    var ident = document.createElement("div");
    ident.className = "ident";
    var name = document.createElement("div");
    name.className = "name mono";
    name.textContent = mutationClass.name;
    ident.appendChild(name);
    var why = document.createElement("div");
    why.className = "why";
    if (mutationClass.kind === "quipu") {
      /* The class name already says "quipu", so spend this line on something
         that differs from row to row: the algebra the class is represented by. */
      why.innerHTML = "";
      var repLabel = document.createElement("span");
      repLabel.textContent = "represented by ";
      why.appendChild(repLabel);
      var rep = document.createElement("code");
      rep.textContent = mutationClass.representative || "A_" + data.length
        + " (no relations)";
      why.appendChild(rep);
      why.title = "the LNA the class is named from, in relation-path notation";
    } else if (mutationClass.kind === "canonical type") {
      why.textContent = "canonical algebra of weight type "
        + mutationClass.form.slice(1) + ", represented by " + mutationClass.representative;
    } else if (mutationClass.kind === "not piecewise hereditary") {
      why.textContent = "in no quipu class at all, represented by "
        + mutationClass.representative;
    } else {
      why.textContent = (mutationClass.form || "no hereditary form reached")
        + ", represented by " + mutationClass.representative;
    }
    ident.appendChild(why);
    summary.appendChild(ident);

    var measure = document.createElement("div");
    measure.className = "measure";
    var pill = document.createElement("span");
    pill.className = "pill " + (PILL[mutationClass.kind] || "other");
    pill.textContent = mutationClass.kind;
    measure.appendChild(pill);
    var bar = document.createElement("div");
    bar.className = "bar";
    var fill = document.createElement("span");
    fill.style.width = Math.max(2, 100 * mutationClass.size / data.classes[0].size) + "%";
    bar.appendChild(fill);
    measure.appendChild(bar);
    var count = document.createElement("div");
    count.className = "count-n";
    count.textContent = mutationClass.size;
    count.title = mutationClass.size + " of " + data.lnaCount + " LNAs";
    measure.appendChild(count);
    summary.appendChild(measure);

    details.appendChild(summary);
    var loaded = false;
    details.addEventListener("toggle", function () {
      if (details.open && !loaded) {
        loaded = true;
        details.appendChild(renderMembers(mutationClass));
      }
    });
    return details;
  }

  function matches(mutationClass, needle) {
    if (kindFilter && mutationClass.kind !== kindFilter) { return false; }
    if (!needle) { return true; }
    if ((mutationClass.name + " " + mutationClass.form + " " + mutationClass.coxeter)
        .toLowerCase().indexOf(needle) >= 0) { return true; }
    return mutationClass.members.some(function (member) {
      return member.relations.indexOf(needle) >= 0;
    });
  }

  function draw() {
    var needle = search.value.trim().toLowerCase();
    var shown = data.classes.filter(function (c) { return matches(c, needle); });
    list.textContent = "";
    if (!shown.length) {
      var empty = document.createElement("p");
      empty.className = "empty";
      empty.textContent = "No class matches. Searching looks at class names, hereditary "
        + "forms, Coxeter polynomials and the relation strings of the members.";
      list.appendChild(empty);
      return;
    }
    shown.forEach(function (c) { list.appendChild(band(c)); });
  }

  chips.forEach(function (chip) {
    chip.addEventListener("click", function () {
      var kind = chip.dataset.kind || null;
      kindFilter = (kindFilter === kind) ? null : kind;
      chips.forEach(function (other) {
        other.setAttribute("aria-pressed", String((other.dataset.kind || null) === kindFilter));
      });
      draw();
    });
  });
  search.addEventListener("input", draw);
  draw();
})();
</script>
"""


def _pluralClasses(count):
    return "1 class" if count == 1 else "{0} classes".format(count)


def render(classification):
    """The page, as the Artifact skeleton wants it: no doctype, head or body."""
    data = collect(classification)
    length = data["length"]
    classCount = len(data["classes"])
    biggest = data["classes"][0] if data["classes"] else None

    chips = ['<button class="chip" aria-pressed="false">all <span class="count">{0}</span></button>'
             .format(classCount)]
    for kind in cv.KINDS:
        counts = data["kinds"].get(kind)
        if not counts:
            continue
        chips.append(
            '<button class="chip" data-kind="{0}" aria-pressed="false">{0} '
            '<span class="count">{1}</span></button>'.format(html.escape(kind), counts[0]))

    collisions = ""
    if data["collisions"]:
        items = []
        for group in data["collisions"]:
            names = "".join("<code>{0}</code>".format(html.escape(n))
                            for n in group["classes"])
            items.append(
                '<div class="collision"><span class="poly mono">{0}</span>'
                '<div class="pair">{1}</div></div>'.format(
                    html.escape(group["polynomial"]), names))
        collisions = """
  <section class="collisions">
    <span class="eyebrow">Where the Coxeter polynomial stops separating</span>
    <h2>{count} shared between classes</h2>
    <p>The Coxeter polynomial is invariant under derived equivalence but does not
       determine it. Each polynomial below is carried by more than one class, so
       nothing about it can tell those classes apart -- they are separated by the
       quipu instead.</p>
    {items}
  </section>""".format(count=("One polynomial" if len(data["collisions"]) == 1
                              else "{0} polynomials".format(len(data["collisions"]))),
                       items="\n    ".join(items))

    unsettled = ""
    if data["unsettled"]:
        unsettled = ('<p class="truncated">{0} LNAs are still unplaced: this '
                     'classification is unfinished.</p>'.format(data["unsettled"]))

    tally = [
        '<span><span class="n">{0}</span> <span class="of">LNAs</span></span>'.format(
            data["lnaCount"]),
        '<span><span class="n">{0}</span> <span class="of">classes</span></span>'.format(
            classCount),
    ]
    if biggest:
        tally.append(
            '<span><span class="n">{0}%</span> <span class="of">of them in the '
            'largest</span></span>'.format(
                round(100 * biggest["size"] / max(data["lnaCount"], 1))))

    return """<title>A_{length} Mutation Classes</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,500;1,6..72,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
{style}
<div class="wrap">
  <header class="page">
    <span class="eyebrow">Linearly oriented Nakayama algebras &middot; derived equivalence</span>
    <h1>A<sub>{length}</sub> up to derived equivalence <span class="sub">&mdash; {classCount}</span></h1>
    <p class="lede">Every linearly oriented Nakayama algebra on {length} vertices, sorted
      into its derived equivalence class. A class named
      <code>P^(m)_(k)</code> is the quipu drawn beside it; open a class to see
      the algebras in it, each as its quiver with an arc over the span of every
      relation, and the mutations that reach it from the representative.</p>
    <div class="tally">{tally}</div>
    {unsettled}
  </header>

  <div class="controls">
    <input type="search" id="search" placeholder="a class, a quipu, a polynomial, or a relation string like 3;4;5"
           aria-label="filter the classes">
    <div class="chips">{chips}</div>
  </div>
{collisions}

  <section class="list" id="list"></section>

  <footer class="page">
    <p>Computed by <code>classify.py {length}</code> and written by
      <code>classes.py {length} --page</code>, from
      <a href="https://arxiv.org/abs/2112.08129">arXiv:2112.08129</a>'s
      combinatorial rule for tilting mutation and the quipu theorem of
      <a href="https://arxiv.org/abs/2305.06642">arXiv:2305.06642</a>.
      Class sizes sum to {lnaCount}, which is Catalan({lengthLess}).</p>
  </footer>
</div>
<script id="data" type="application/json">{json}</script>
<script>window.__classification = JSON.parse(document.getElementById("data").textContent);</script>
{script}
""".format(
        length=length,
        lengthLess=length - 1,
        classCount=_pluralClasses(classCount),
        lnaCount=data["lnaCount"],
        style=STYLE,
        tally="".join(tally),
        chips="".join(chips),
        collisions=collisions,
        unsettled=unsettled,
        json=json.dumps(data, separators=(",", ":")).replace("</", "<\\/"),
        script=SCRIPT,
    )


def renderStandalone(classification):
    """The same page as a file to open in a browser."""
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '</head>\n<body>\n'
        + render(classification)
        + '</body>\n</html>\n'
    )
