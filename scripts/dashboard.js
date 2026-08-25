/* GDSC RNA-Seq dashboard runtime.
   No dependencies: every chart is hand-built SVG so the report stays offline-safe.
   Sample colours come from CSS custom properties, so the dark-mode toggle repaints
   the charts without a re-render. */
(function () {
"use strict";

var D = JSON.parse(document.getElementById("payload").textContent);
var SAMPLES = D.samples || [];
var M = D.metrics || {};
var NPAL = 12;

/* ---------------------------------------------------------------- helpers */
function $(sel, root) { return (root || document).querySelector(sel); }
function el(tag, attrs, kids) {
  var n = document.createElement(tag);
  for (var k in attrs || {}) {
    if (k === "class") n.className = attrs[k];
    else if (k === "html") n.innerHTML = attrs[k];
    else if (k === "text") n.textContent = attrs[k];
    else n.setAttribute(k, attrs[k]);
  }
  (kids || []).forEach(function (c) { n.appendChild(c); });
  return n;
}
function svg(tag, attrs) {
  var n = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (var k in attrs || {}) {
    if (k === "text") n.textContent = attrs[k];
    else n.setAttribute(k, attrs[k]);
  }
  return n;
}
function color(i) { return "var(--c" + ((i % NPAL) + 1) + ")"; }
function sampleColor(s) { return color(Math.max(0, SAMPLES.indexOf(s))); }

// `groups` is present only when the sample sheet carried a group column. When it is,
// the PCA is coloured by group rather than by sample, so replicates read as one set.
var GROUPS = (D.pca && D.pca.groups) || null;
var GROUP_LEVELS = GROUPS
  ? Array.from(new Set(GROUPS.filter(Boolean))).sort(function (a, b) {
      return a.localeCompare(b, undefined, { numeric: true });
    })
  : [];
function groupOf(name) {
  if (!GROUPS || !D.pca) return null;
  var i = D.pca.samples.indexOf(name);
  return i < 0 ? null : GROUPS[i];
}
function groupColor(g) {
  var i = GROUP_LEVELS.indexOf(g);
  return i < 0 ? "var(--text-3)" : color(i);
}
function pcaColor(name) {
  if (!GROUPS) return sampleColor(name);
  var g = groupOf(name);
  return g ? groupColor(g) : "var(--text-3)";
}

function fmt(v, kind) {
  if (v === null || v === undefined || (typeof v === "number" && !isFinite(v))) return null;
  switch (kind) {
    case "reads":
      if (v >= 1e9) return (v / 1e9).toFixed(2) + "B";
      if (v >= 1e6) return (v / 1e6).toFixed(2) + "M";
      if (v >= 1e3) return (v / 1e3).toFixed(1) + "K";
      return String(Math.round(v));
    case "pct":  return v.toFixed(1) + "%";
    case "pct1": return v.toFixed(1) + "%";
    case "num2": return v.toFixed(2);
    default:
      if (Math.abs(v) >= 1000) return v.toLocaleString(undefined, { maximumFractionDigits: 0 });
      return String(Math.round(v * 100) / 100);
  }
}
function median(a) {
  var b = a.filter(function (x) { return x !== null && isFinite(x); }).sort(function (x, y) { return x - y; });
  if (!b.length) return null;
  var m = Math.floor(b.length / 2);
  return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2;
}
function niceTicks(lo, hi, want) {
  if (hi === lo) { hi = lo + 1; }
  var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / (want || 5))));
  var err = span / (want || 5) / step;
  if (err >= 7.5) step *= 10; else if (err >= 3.5) step *= 5; else if (err >= 1.5) step *= 2;
  var out = [], t = Math.ceil(lo / step) * step;
  for (; t <= hi + step * 1e-9; t += step) out.push(Math.round(t / step) * step);
  return out;
}

/* --------------------------------------------------------------- tooltip */
var tipEl = $("#tip");
function showTip(evt, title, rows) {
  var h = "<b>" + esc(title) + "</b>";
  (rows || []).forEach(function (r) {
    h += '<div class="row"><span>' + esc(r[0]) + "</span><span>" + esc(r[1]) + "</span></div>";
  });
  tipEl.innerHTML = h;
  tipEl.style.opacity = "1";
  moveTip(evt);
}
function moveTip(evt) {
  var pad = 14, w = tipEl.offsetWidth, h = tipEl.offsetHeight;
  var x = evt.clientX + pad, y = evt.clientY + pad;
  if (x + w > innerWidth - 8) x = evt.clientX - w - pad;
  if (y + h > innerHeight - 8) y = evt.clientY - h - pad;
  tipEl.style.left = Math.max(8, x) + "px";
  tipEl.style.top = Math.max(8, y) + "px";
}
function hideTip() { tipEl.style.opacity = "0"; }
function esc(s) {
  return String(s).replace(/[&<>"]/g, function (c) {
    return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
  });
}

/* ------------------------------------------------- cross-panel highlight */
var hovered = null, pinned = null, groupHover = null;
function current() { return hovered || pinned; }
function setHover(s) { hovered = s; applyHighlight(); }
function setGroupHover(list) { groupHover = list; applyHighlight(); }
function togglePin(s) { pinned = (pinned === s) ? null : s; applyHighlight(); }
function applyHighlight() {
  var cur = current();
  // A hovered group highlights every one of its samples; a hovered sample wins over it.
  var set = cur ? [cur] : groupHover;
  document.querySelectorAll("[data-sample]").forEach(function (n) {
    var mine = set && set.indexOf(n.getAttribute("data-sample")) >= 0;
    n.classList.toggle("active", !!set && mine);
    n.classList.toggle("dim", !!set && !mine);
  });
}
function bindSample(node, s, title, rows) {
  node.setAttribute("data-sample", s);
  node.addEventListener("mouseenter", function (e) { setHover(s); showTip(e, title, rows()); });
  node.addEventListener("mousemove", moveTip);
  node.addEventListener("mouseleave", function () { setHover(null); hideTip(); });
  node.addEventListener("click", function () { togglePin(s); });
  node.style.cursor = "pointer";
}

/* ============================================================== BAR PANEL */
var barSort = "sample";
var barPanels = [];

function barPanel(def) {
  var card = el("div", { class: "card" });
  var head = el("div", { class: "card-head" });
  head.appendChild(el("div", {}, [
    el("h3", { text: def.label }),
    el("div", { class: "hint", text: def.desc || "" })
  ]));
  card.appendChild(head);
  var host = el("div", { class: "card-pad" });
  card.appendChild(host);

  function draw() {
    host.innerHTML = "";
    var rows = SAMPLES.map(function (s) { return { s: s, v: M[s] ? M[s][def.key] : null }; })
      .filter(function (r) { return r.v !== null && r.v !== undefined && isFinite(r.v); });
    if (!rows.length) { host.appendChild(el("div", { class: "empty", text: "No data" })); return; }
    if (barSort === "value") rows.sort(function (a, b) { return b.v - a.v; });

    var W = 480, H = 262, mL = 52, mR = 10, mT = 12, mB = 60;
    var iw = W - mL - mR, ih = H - mT - mB;
    var isPct = def.fmt === "pct" || def.fmt === "pct1";
    var maxV = Math.max.apply(null, rows.map(function (r) { return r.v; }));
    var hi = isPct ? Math.min(100, Math.max(10, Math.ceil(maxV / 10) * 10)) : maxV * 1.08;
    if (hi <= 0) hi = 1;
    var med = median(rows.map(function (r) { return r.v; }));

    var s = svg("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, role: "img",
      "aria-label": def.label + " by sample" });
    var y = function (v) { return mT + ih - (v / hi) * ih; };

    var g = svg("g", { class: "grid" });
    niceTicks(0, hi, 5).forEach(function (t) {
      g.appendChild(svg("line", { x1: mL, x2: mL + iw, y1: y(t), y2: y(t) }));
      g.appendChild(svg("text", { x: mL - 7, y: y(t) + 3.5, "text-anchor": "end",
        text: def.fmt === "reads" ? fmt(t, "reads") : (Math.round(t * 100) / 100) }));
    });
    s.appendChild(g);

    var bw = iw / rows.length, pad = Math.min(5, bw * 0.22);
    rows.forEach(function (r, i) {
      var x = mL + i * bw + pad / 2, w = Math.max(1, bw - pad);
      // Bars stay a single accent colour: sample identity is carried by the
      // highlight/dim link, and 20 hues here would read as noise.
      var rect = svg("rect", { class: "barrect", x: x, y: y(r.v), width: w,
        height: Math.max(1, mT + ih - y(r.v)), rx: 2 });
      bindSample(rect, r.s, r.s, function () {
        return [[def.label, fmt(r.v, def.fmt)], ["Cohort median", fmt(med, def.fmt)]];
      });
      s.appendChild(rect);
    });

    if (med !== null) {
      s.appendChild(svg("line", { class: "refline", x1: mL, x2: mL + iw, y1: y(med), y2: y(med) }));
    }

    var ax = svg("g", { class: "axis" });
    ax.appendChild(svg("line", { x1: mL, x2: mL + iw, y1: mT + ih, y2: mT + ih }));
    s.appendChild(ax);

    // Only label the x axis while the names still fit; otherwise the count line does.
    if (rows.length <= 22) {
      rows.forEach(function (r, i) {
        var cx = mL + i * bw + bw / 2;
        var t = svg("text", { x: cx, y: mT + ih + 6, "text-anchor": "end",
          transform: "rotate(-52 " + cx + " " + (mT + ih + 6) + ")", text: r.s });
        t.setAttribute("data-sample", r.s);
        s.appendChild(t);
      });
    } else {
      s.appendChild(svg("text", { x: mL + iw / 2, y: H - 12, "text-anchor": "middle",
        text: rows.length + " samples" }));
    }
    host.appendChild(s);
  }

  barPanels.push(draw);
  draw();
  return card;
}

/* ========================================================= COVERAGE LINES */
function coveragePanel(host) {
  var cov = D.coverage || {};
  var keys = Object.keys(cov).filter(function (k) { return SAMPLES.indexOf(k) >= 0; });
  host.innerHTML = "";
  host.appendChild(el("div", { class: "card-head" }, [
    el("div", {}, [
      el("h3", { text: "Normalised coverage along the transcript body" }),
      el("div", { class: "hint", text: "Picard CollectRnaSeqMetrics histogram · " + keys.length + " samples" })
    ])
  ]));
  if (!keys.length) {
    host.appendChild(el("div", { class: "empty", text: "No Picard coverage histogram found." }));
    return;
  }

  var W = 960, H = 320, mL = 52, mR = 14, mT = 14, mB = 38;
  var iw = W - mL - mR, ih = H - mT - mB;
  var maxY = 0, maxX = 0;
  keys.forEach(function (k) {
    cov[k].forEach(function (p) {
      if (p[1] > maxY) maxY = p[1];
      if (p[0] > maxX) maxX = p[0];
    });
  });
  maxY = Math.max(0.5, maxY * 1.08);
  // Picard reports normalized_position as 0-100; older/other emitters use 0-1.
  var xspan = maxX > 1.5 ? maxX : 1;

  var x = function (v) { return mL + (v / xspan) * iw; };
  var y = function (v) { return mT + ih - (v / maxY) * ih; };
  var s = svg("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, role: "img",
    "aria-label": "Gene body coverage" });

  var g = svg("g", { class: "grid" });
  niceTicks(0, maxY, 5).forEach(function (t) {
    g.appendChild(svg("line", { x1: mL, x2: mL + iw, y1: y(t), y2: y(t) }));
    g.appendChild(svg("text", { x: mL - 7, y: y(t) + 3.5, "text-anchor": "end",
      text: Math.round(t * 100) / 100 }));
  });
  [0, 0.25, 0.5, 0.75, 1].forEach(function (t) {
    g.appendChild(svg("line", { x1: x(t * xspan), x2: x(t * xspan), y1: mT, y2: mT + ih }));
    g.appendChild(svg("text", { x: x(t * xspan), y: H - 20, "text-anchor": "middle",
      text: t.toFixed(2) }));
  });
  s.appendChild(g);
  s.appendChild(svg("text", { x: mL, y: H - 6, "text-anchor": "start", text: "5' end" }));
  s.appendChild(svg("text", { x: mL + iw, y: H - 6, "text-anchor": "end", text: "3' end" }));
  s.appendChild(svg("text", { x: 12, y: mT + ih / 2, "text-anchor": "middle",
    transform: "rotate(-90 12 " + (mT + ih / 2) + ")", text: "Normalised coverage" }));

  keys.forEach(function (k) {
    var pts = cov[k].map(function (p) { return x(p[0]) .toFixed(1) + "," + y(p[1]).toFixed(1); }).join(" ");
    var line = svg("polyline", { class: "cov", points: pts, stroke: sampleColor(k) });
    var peak = cov[k].reduce(function (a, b) { return b[1] > a[1] ? b : a; }, cov[k][0]);
    bindSample(line, k, k, function () {
      return [["Peak coverage", (Math.round(peak[1] * 100) / 100)],
              ["at position", (peak[0] / xspan).toFixed(2)]];
    });
    s.appendChild(line);
  });
  host.appendChild(el("div", { class: "card-pad" }, [s]));
  host.appendChild(sampleLegend(keys));
}

function sampleLegend(keys) {
  var wrap = el("div", { class: "legend" });
  keys.forEach(function (k) {
    var sw = el("i");
    sw.style.background = sampleColor(k);
    var sp = el("span", {}, [sw, el("span", { text: k })]);
    sp.setAttribute("data-sample", k);
    sp.addEventListener("mouseenter", function () { setHover(k); });
    sp.addEventListener("mouseleave", function () { setHover(null); });
    sp.addEventListener("click", function () { togglePin(k); });
    wrap.appendChild(sp);
  });
  return wrap;
}

/* ================================================================== PCA */
var pcaX = 1, pcaY = 2;

function pcaPanel(host) {
  var p = D.pca;
  host.innerHTML = "";
  host.appendChild(el("div", { class: "card-head" }, [
    el("div", {}, [
      el("h3", { text: "Principal component analysis" }),
      el("div", { class: "hint", text: p ? (p.n_hvg.toLocaleString() + " variable genes of " +
        p.n_genes_total.toLocaleString() + " · " + p.hvg_method) : "" })
    ])
  ]));
  if (!p) {
    host.appendChild(el("div", { class: "empty",
      text: (D.sources && D.sources.counts) ? "PCA needs at least 2 samples and 2 genes."
                                            : "No PCA outputs or counts matrix found." }));
    return;
  }

  var W = 620, H = 470, mL = 58, mR = 18, mT = 16, mB = 50;
  var iw = W - mL - mR, ih = H - mT - mB;
  var xi = pcaX - 1, yi = pcaY - 1;
  var xs = p.scores.map(function (r) { return r[xi]; });
  var ys = p.scores.map(function (r) { return r[yi]; });
  var xr = pad(Math.min.apply(null, xs), Math.max.apply(null, xs));
  var yr = pad(Math.min.apply(null, ys), Math.max.apply(null, ys));
  function pad(lo, hi) { var d = (hi - lo) || 1; return [lo - d * 0.12, hi + d * 0.12]; }
  var X = function (v) { return mL + (v - xr[0]) / (xr[1] - xr[0]) * iw; };
  var Y = function (v) { return mT + ih - (v - yr[0]) / (yr[1] - yr[0]) * ih; };

  var s = svg("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, role: "img",
    "aria-label": "PCA scatter plot" });
  var g = svg("g", { class: "grid" });
  niceTicks(xr[0], xr[1], 5).forEach(function (t) {
    g.appendChild(svg("line", { x1: X(t), x2: X(t), y1: mT, y2: mT + ih }));
    g.appendChild(svg("text", { x: X(t), y: mT + ih + 15, "text-anchor": "middle", text: t }));
  });
  niceTicks(yr[0], yr[1], 5).forEach(function (t) {
    g.appendChild(svg("line", { x1: mL, x2: mL + iw, y1: Y(t), y2: Y(t) }));
    g.appendChild(svg("text", { x: mL - 8, y: Y(t) + 3.5, "text-anchor": "end", text: t }));
  });
  s.appendChild(g);
  if (xr[0] < 0 && xr[1] > 0) s.appendChild(svg("line", { class: "refline", x1: X(0), x2: X(0), y1: mT, y2: mT + ih }));
  if (yr[0] < 0 && yr[1] > 0) s.appendChild(svg("line", { class: "refline", x1: mL, x2: mL + iw, y1: Y(0), y2: Y(0) }));

  var showLabels = p.samples.length <= 30;
  p.samples.forEach(function (name, i) {
    var cx = X(p.scores[i][xi]), cy = Y(p.scores[i][yi]);
    if (showLabels) {
      var t = svg("text", { class: "ptext", x: cx + 9, y: cy - 7, text: name });
      t.setAttribute("data-sample", name);
      s.appendChild(t);
    }
    var c = svg("circle", { class: "pt", cx: cx, cy: cy, r: 6.5, fill: pcaColor(name) });
    bindSample(c, name, name, function () {
      var rows = [];
      if (GROUPS && GROUPS[i]) rows.push(["Group", GROUPS[i]]);
      rows.push(["PC" + pcaX, p.scores[i][xi].toFixed(2)]);
      rows.push(["PC" + pcaY, p.scores[i][yi].toFixed(2)]);
      return rows;
    });
    s.appendChild(c);
  });

  s.appendChild(svg("text", { x: mL + iw / 2, y: H - 12, "text-anchor": "middle",
    text: "PC" + pcaX + " — " + p.var_explained[xi].toFixed(1) + "% of variance" }));
  s.appendChild(svg("text", { x: 14, y: mT + ih / 2, "text-anchor": "middle",
    transform: "rotate(-90 14 " + (mT + ih / 2) + ")",
    text: "PC" + pcaY + " — " + p.var_explained[yi].toFixed(1) + "% of variance" }));

  host.appendChild(el("div", { class: "card-pad" }, [s]));
  if (GROUP_LEVELS.length) host.appendChild(groupLegend());
}

function groupLegend() {
  var wrap = el("div", { class: "legend" });
  GROUP_LEVELS.forEach(function (g) {
    var sw = el("i");
    sw.style.background = groupColor(g);
    var members = D.pca.samples.filter(function (nm) { return groupOf(nm) === g; });
    var sp = el("span", { title: members.join(", ") }, [
      sw, el("span", { text: g + " (" + members.length + ")" })
    ]);
    // A group is many samples, so drive the shared highlight for each of them.
    sp.addEventListener("mouseenter", function () { setGroupHover(members); });
    sp.addEventListener("mouseleave", function () { setGroupHover(null); });
    wrap.appendChild(sp);
  });
  if (GROUPS.some(function (g) { return !g; })) {
    var sw = el("i");
    sw.style.background = "var(--text-3)";
    wrap.appendChild(el("span", { title: "no group in the sample sheet" },
      [sw, el("span", { text: "ungrouped" })]));
  }
  return wrap;
}

function screePanel(host) {
  var p = D.pca;
  host.innerHTML = "";
  host.appendChild(el("div", { class: "card-head" }, [
    el("div", {}, [
      el("h3", { text: "Variance explained" }),
      el("div", { class: "hint", text: "Per component, with cumulative total" })
    ])
  ]));
  if (!p) { host.appendChild(el("div", { class: "empty", text: "—" })); return; }

  var v = p.var_explained, W = 460, H = 470, mL = 46, mR = 40, mT = 16, mB = 50;
  var iw = W - mL - mR, ih = H - mT - mB;
  var hi = Math.max.apply(null, v) * 1.12;
  var y = function (t) { return mT + ih - t / hi * ih; };
  var yc = function (t) { return mT + ih - t / 100 * ih; };

  var s = svg("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, role: "img",
    "aria-label": "Scree plot" });
  var g = svg("g", { class: "grid" });
  niceTicks(0, hi, 5).forEach(function (t) {
    g.appendChild(svg("line", { x1: mL, x2: mL + iw, y1: y(t), y2: y(t) }));
    g.appendChild(svg("text", { x: mL - 7, y: y(t) + 3.5, "text-anchor": "end", text: t + "%" }));
  });
  s.appendChild(g);

  var bw = iw / v.length, cum = 0, cumPts = [];
  v.forEach(function (val, i) {
    var x = mL + i * bw + bw * 0.16, w = bw * 0.68;
    var rect = svg("rect", { class: "barrect", x: x, y: y(val), width: w,
      height: Math.max(1, mT + ih - y(val)), rx: 2, fill: (i + 1 === pcaX || i + 1 === pcaY) ? "var(--accent)" : "var(--text-3)" });
    rect.style.cursor = "pointer";
    rect.addEventListener("mouseenter", function (e) {
      showTip(e, "PC" + (i + 1), [["Variance", val.toFixed(2) + "%"]]);
    });
    rect.addEventListener("mousemove", moveTip);
    rect.addEventListener("mouseleave", hideTip);
    rect.addEventListener("click", function () { pcaY = pcaX; pcaX = i + 1; renderPCA(); });
    s.appendChild(rect);
    cum += val;
    cumPts.push((mL + i * bw + bw / 2).toFixed(1) + "," + yc(cum).toFixed(1));
    s.appendChild(svg("text", { x: mL + i * bw + bw / 2, y: mT + ih + 15, "text-anchor": "middle",
      text: "PC" + (i + 1) }));
  });
  s.appendChild(svg("polyline", { points: cumPts.join(" "), fill: "none",
    stroke: "var(--c3)", "stroke-width": 1.6, "stroke-dasharray": "4 3" }));
  s.appendChild(svg("text", { x: mL + iw, y: yc(cum) - 8, "text-anchor": "end",
    fill: "var(--c3)", text: "cumulative " + cum.toFixed(0) + "%" }));
  s.appendChild(svg("text", { x: mL + iw / 2, y: H - 10, "text-anchor": "middle",
    text: "Click a bar to plot that component" }));
  host.appendChild(el("div", { class: "card-pad" }, [s]));
}

function renderPCA() {
  pcaPanel($("#pcaCard"));
  screePanel($("#screeCard"));
  applyHighlight();
  buildPcaControls();
}

function buildPcaControls() {
  var host = $("#pcaCtl");
  host.innerHTML = "";
  if (!D.pca) return;
  var n = D.pca.n_comp;
  function sel(id, val, onchange) {
    var s = el("select", { id: id, "aria-label": id });
    for (var i = 1; i <= n; i++) {
      var o = el("option", { value: i, text: "PC" + i });
      if (i === val) o.selected = true;
      s.appendChild(o);
    }
    s.addEventListener("change", function () { onchange(parseInt(s.value, 10)); });
    return s;
  }
  host.appendChild(el("span", { class: "chip", text: "x-axis" }));
  host.appendChild(sel("pcx", pcaX, function (v) { pcaX = v; renderPCA(); }));
  host.appendChild(el("span", { class: "chip", text: "y-axis" }));
  host.appendChild(sel("pcy", pcaY, function (v) { pcaY = v; renderPCA(); }));
}

/* ==================================================== CORRELATION HEATMAP */
function parseColor(str) {
  var s = String(str).trim();
  var m = s.match(/^#([0-9a-f]{3}|[0-9a-f]{6})$/i);
  if (m) {
    var h = m[1];
    if (h.length === 3) h = h[0] + h[0] + h[1] + h[1] + h[2] + h[2];
    return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
  }
  m = s.match(/rgba?\(([^)]+)\)/i);
  if (m) return m[1].split(",").slice(0, 3).map(function (v) { return parseInt(v, 10) || 0; });
  return [128, 128, 128];
}

// Read the ramp from the theme tokens so the heatmap follows dark mode instead of
// burning a light-mode palette into the cells. Distances are sequential, not
// diverging: 0 is a real floor (a sample against itself), so a two-ended red/blue
// scale would imply a midpoint that does not exist.
function heatStops() {
  var cs = getComputedStyle(document.documentElement);
  return [parseColor(cs.getPropertyValue("--dist-near")),
          parseColor(cs.getPropertyValue("--dist-far"))];
}

function heatColor(t, stops) {
  var st = stops || heatStops();
  var a = st[0], b = st[1];
  return "rgb(" + a.map(function (c, i) { return Math.round(c + (b[i] - c) * t); }).join(",") + ")";
}

function distPanel(host) {
  var c = D.dist;
  host.innerHTML = "";
  host.appendChild(el("div", { class: "card-head" }, [
    el("div", {}, [
      el("h3", { text: "Sample-to-sample distance" }),
      el("div", { class: "hint", text: c ? (c.label + " on the normalised matrix used for "
        + "the PCA (" + c.n_genes.toLocaleString() + " genes)") : "" })
    ])
  ]));
  if (!c) {
    host.appendChild(el("div", { class: "empty",
      text: "No normalised count matrix — sample distances omitted." }));
    return;
  }

  var n = c.samples.length;
  var off = [];
  for (var i = 0; i < n; i++) for (var j = 0; j < n; j++) if (i !== j) off.push(c.matrix[i][j]);
  // Stretch the ramp across the off-diagonal range. Anchoring at 0 would spend most of
  // the scale on the empty gap below the closest pair and flatten the real structure;
  // the diagonal is painted at the near endpoint regardless.
  var lo = Math.min.apply(null, off), hi = Math.max.apply(null, off);
  if (hi - lo < 1e-9) { lo = 0; hi = hi || 1; }
  var norm = function (v) { return Math.max(0, Math.min(1, (v - lo) / (hi - lo))); };

  var stops = heatStops();
  var labels = n <= 40;
  var pad = labels ? Math.min(120, 8 + Math.max.apply(null, c.samples.map(function (s) { return s.length; })) * 5.4) : 10;
  var cell = Math.max(6, Math.min(30, 430 / n));
  var W = pad + cell * n + 12, H = pad + cell * n + 12;

  var s = svg("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, role: "img",
    "aria-label": "Sample distance heatmap" });
  c.samples.forEach(function (rs, ri) {
    c.samples.forEach(function (cs, ci) {
      var v = c.matrix[ri][ci];
      var r = svg("rect", { x: pad + ci * cell, y: pad + ri * cell,
        width: cell - 0.6, height: cell - 0.6, rx: 1.5,
        fill: ri === ci ? heatColor(0, stops) : heatColor(norm(v), stops) });
      r.style.cursor = "pointer";
      r.addEventListener("mouseenter", function (e) {
        setHover(ri === ci ? rs : null);
        showTip(e, rs + "  \u00d7  " + cs, [[c.label, v.toFixed(3)]]);
      });
      r.addEventListener("mousemove", moveTip);
      r.addEventListener("mouseleave", function () { setHover(null); hideTip(); });
      s.appendChild(r);
    });
    if (labels) {
      var ty = svg("text", { x: pad - 6, y: pad + ri * cell + cell / 2 + 3,
        "text-anchor": "end", text: rs });
      ty.setAttribute("data-sample", rs);
      s.appendChild(ty);
      var tx = svg("text", { x: pad + ri * cell + cell / 2, y: pad - 6, "text-anchor": "start",
        transform: "rotate(-90 " + (pad + ri * cell + cell / 2) + " " + (pad - 6) + ")", text: rs });
      tx.setAttribute("data-sample", rs);
      s.appendChild(tx);
    }
  });

  var body = el("div", { class: "card-pad" }, [s]);
  var ramp = el("div", { class: "ramp" });
  ramp.style.background = "linear-gradient(90deg," + heatColor(0, stops) + "," +
    heatColor(0.5, stops) + "," + heatColor(1, stops) + ")";
  body.appendChild(el("div", { class: "scale", style: "margin-top:10px" }, [
    el("span", { text: lo.toFixed(1) + " · closest" }), ramp,
    el("span", { text: hi.toFixed(1) + " · most distant" })
  ]));
  host.appendChild(body);
}

/* ================================================== GENE VARIANCE CURVE */
function varPanel(host) {
  var p = D.pca;
  host.innerHTML = "";
  host.appendChild(el("div", { class: "card-head" }, [
    el("div", {}, [
      el("h3", { text: "Gene variance by rank" }),
      el("div", { class: "hint", text: p ? ("dashed line = " + p.n_hvg.toLocaleString() +
        " genes used for PCA") : "" })
    ])
  ]));
  if (!p || !p.var_curve || !p.var_curve.length) {
    host.appendChild(el("div", { class: "empty", text: "—" })); return;
  }

  var pts = p.var_curve, W = 620, H = 300, mL = 52, mR = 16, mT = 16, mB = 42;
  var iw = W - mL - mR, ih = H - mT - mB;
  var maxR = pts[pts.length - 1][0], maxV = Math.max.apply(null, pts.map(function (q) { return q[1]; })) * 1.06;
  var X = function (r) { return mL + Math.log10(Math.max(1, r)) / Math.log10(maxR) * iw; };
  var Y = function (v) { return mT + ih - v / (maxV || 1) * ih; };

  var s = svg("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, role: "img",
    "aria-label": "Gene variance versus rank" });
  var g = svg("g", { class: "grid" });
  niceTicks(0, maxV, 4).forEach(function (t) {
    g.appendChild(svg("line", { x1: mL, x2: mL + iw, y1: Y(t), y2: Y(t) }));
    g.appendChild(svg("text", { x: mL - 7, y: Y(t) + 3.5, "text-anchor": "end",
      text: Math.round(t * 100) / 100 }));
  });
  var dec = Math.ceil(Math.log10(maxR));
  for (var e = 0; e <= dec; e++) {
    var r = Math.pow(10, e);
    if (r > maxR) break;
    g.appendChild(svg("line", { x1: X(r), x2: X(r), y1: mT, y2: mT + ih }));
    g.appendChild(svg("text", { x: X(r), y: mT + ih + 15, "text-anchor": "middle",
      text: r >= 1000 ? (r / 1000) + "k" : r }));
  }
  s.appendChild(g);
  s.appendChild(svg("polyline", {
    points: pts.map(function (q) { return X(q[0]).toFixed(1) + "," + Y(q[1]).toFixed(1); }).join(" "),
    fill: "none", stroke: "var(--accent)", "stroke-width": 1.8
  }));
  s.appendChild(svg("line", { class: "refline", x1: X(p.n_hvg), x2: X(p.n_hvg), y1: mT, y2: mT + ih }));
  s.appendChild(svg("text", { x: mL + iw / 2, y: H - 8, "text-anchor": "middle",
    text: "Gene rank (log scale)" }));
  host.appendChild(el("div", { class: "card-pad" }, [s]));
}

/* ==================================================== SAMPLE METRIC TABLE */
var metricSort = { key: "sample", dir: 1 };

function activeDefs() {
  var pick = $("#groupFilter").value;
  return (D.metric_defs || []).filter(function (d) {
    if (pick !== "__all__" && d.group !== pick) return false;
    return SAMPLES.some(function (s) {
      var v = M[s] ? M[s][d.key] : null;
      return v !== null && v !== undefined && isFinite(v);
    });
  });
}

function renderMetricTable() {
  var defs = activeDefs();
  var tbl = $("#metricTable");
  tbl.innerHTML = "";
  if (!defs.length) {
    tbl.appendChild(el("tbody", {}, [el("tr", {}, [el("td", { class: "empty", text: "No metrics parsed." })])]));
    return;
  }

  // Column maxima drive the in-cell proportion bars.
  var maxes = {};
  defs.forEach(function (d) {
    var vals = SAMPLES.map(function (s) { return M[s] ? M[s][d.key] : null; })
      .filter(function (v) { return v !== null && v !== undefined && isFinite(v); });
    maxes[d.key] = vals.length ? Math.max.apply(null, vals) : 0;
  });

  var thead = el("thead");
  var groupRow = el("tr", { class: "groups" });
  groupRow.appendChild(el("th", { class: "sticky-l", text: "" }));
  var run = null;
  defs.forEach(function (d) {
    if (run && run.g === d.group) { run.n++; return; }
    run = { g: d.group, n: 1, th: el("th", { text: d.group }) };
    groupRow.appendChild(run.th);
    run.th.__run = run;
  });
  Array.prototype.forEach.call(groupRow.querySelectorAll("th"), function (th) {
    if (th.__run) th.setAttribute("colspan", th.__run.n);
  });
  thead.appendChild(groupRow);

  var colRow = el("tr", { class: "cols" });
  var th0 = el("th", { class: "sticky-l" , html: "Sample<span class='arrow'>▼</span>" });
  th0.addEventListener("click", function () { sortMetrics("sample"); });
  if (metricSort.key === "sample") th0.classList.add("sorted");
  colRow.appendChild(th0);
  defs.forEach(function (d) {
    var th = el("th", { title: d.desc || d.label, html: esc(d.label) + "<span class='arrow'>▼</span>" });
    if (metricSort.key === d.key) th.classList.add("sorted");
    th.addEventListener("click", function () { sortMetrics(d.key); });
    colRow.appendChild(th);
  });
  thead.appendChild(colRow);
  tbl.appendChild(thead);

  var q = $("#metricSearch").value.trim().toLowerCase();
  var rows = SAMPLES.slice();
  rows.sort(function (a, b) {
    if (metricSort.key === "sample") return a.localeCompare(b) * metricSort.dir;
    var va = M[a] ? M[a][metricSort.key] : null, vb = M[b] ? M[b][metricSort.key] : null;
    if (va === null || va === undefined) return 1;
    if (vb === null || vb === undefined) return -1;
    return (va - vb) * metricSort.dir;
  });

  var tbody = el("tbody");
  rows.forEach(function (s) {
    if (q && s.toLowerCase().indexOf(q) < 0) return;
    var tr = el("tr");
    tr.setAttribute("data-sample", s);
    var td0 = el("td", { class: "sticky-l" });
    var sw = el("i");
    sw.style.cssText = "display:inline-block;width:8px;height:8px;border-radius:2px;margin-right:8px;background:" + pcaColor(s);
    td0.appendChild(sw);
    td0.appendChild(document.createTextNode(s));
    tr.appendChild(td0);
    defs.forEach(function (d) {
      var v = M[s] ? M[s][d.key] : null;
      var td = el("td");
      var txt = fmt(v, d.fmt);
      if (txt === null) { td.className = "na"; td.textContent = "—"; }
      else {
        var span = el("span", { class: "cell", text: txt });
        if (maxes[d.key] > 0 && v > 0) {
          var bar = el("span", { class: "bar" });
          bar.style.transform = "scaleX(" + Math.min(1, v / maxes[d.key]).toFixed(3) + ")";
          span.appendChild(bar);
        }
        td.appendChild(span);
      }
      tr.appendChild(td);
    });
    tr.addEventListener("mouseenter", function () { setHover(s); });
    tr.addEventListener("mouseleave", function () { setHover(null); });
    tr.addEventListener("click", function () { togglePin(s); });
    tbody.appendChild(tr);
  });
  tbl.appendChild(tbody);
  applyHighlight();
}

function sortMetrics(key) {
  if (metricSort.key === key) metricSort.dir *= -1;
  else metricSort = { key: key, dir: key === "sample" ? 1 : -1 };
  renderMetricTable();
}

/* ========================================================== GENE TABLE */
var geneUnit = "tpm", geneMito = "all", geneSort = { key: "mean", dir: -1 };

function geneValues(row) { return geneUnit === "counts" ? (row.counts || []) : row.tpm; }
function geneMean(row) { return geneUnit === "counts" ? (row.mean_count || 0) : row.mean_tpm; }

function renderGeneTable() {
  var g = D.top_genes;
  var tbl = $("#geneTable");
  tbl.innerHTML = "";
  if (!g || !g.rows.length) {
    tbl.appendChild(el("tbody", {}, [el("tr", {}, [el("td", { class: "empty",
      text: "No TPM matrix found (featurecounts/featurecounts.readcounts_tpm.ann.tsv)." })])]));
    return;
  }
  if (geneUnit === "counts" && !g.has_counts) geneUnit = "tpm";

  var samples = g.samples;
  var thead = el("thead");
  var tr = el("tr", { class: "cols" });
  function hdr(label, key, cls, title) {
    // key === null marks a fixed column: no arrow, no click, no pointer cursor.
    var sortable = key !== null;
    var th = el("th", {
      class: (cls || "") + (sortable ? "" : " nosort"),
      title: title || label,
      html: esc(label) + (sortable ? "<span class='arrow'>▼</span>" : "")
    });
    if (!sortable) return th;
    if (geneSort.key === key) th.classList.add("sorted");
    th.addEventListener("click", function () {
      if (geneSort.key === key) geneSort.dir *= -1;
      else geneSort = { key: key, dir: -1 };
      renderGeneTable();
    });
    return th;
  }
  tr.appendChild(hdr("Gene", null, "sticky-l", "Search the box above to find a gene"));
  tr.appendChild(hdr("Chr", "chr"));
  tr.appendChild(hdr(geneUnit === "counts" ? "Mean count" : "Mean TPM", "mean"));
  samples.forEach(function (s, i) {
    var th = hdr(s, "s" + i);
    th.setAttribute("data-sample", s);
    tr.appendChild(th);
  });
  thead.appendChild(tr);
  tbl.appendChild(thead);

  var q = $("#geneSearch").value.trim().toLowerCase();
  var rows = g.rows.filter(function (r) {
    if (geneMito === "mito" && !r.mito) return false;
    if (geneMito === "nomito" && r.mito) return false;
    if (q && r.gene.toLowerCase().indexOf(q) < 0 && r.id.toLowerCase().indexOf(q) < 0) return false;
    return true;
  });
  rows.sort(function (a, b) {
    var k = geneSort.key;
    if (k === "chr") return String(a.chr).localeCompare(String(b.chr), undefined, { numeric: true }) * geneSort.dir;
    if (k === "mean") return (geneMean(a) - geneMean(b)) * geneSort.dir;
    var i = parseInt(k.slice(1), 10);
    return ((geneValues(a)[i] || 0) - (geneValues(b)[i] || 0)) * geneSort.dir;
  });

  // Shade each cell against its own gene's maximum. Scaling to the table-wide
  // maximum would wash every row out, since the top genes span three orders of
  // magnitude; row-relative shading is what makes across-sample differences visible.

  var tbody = el("tbody");
  rows.forEach(function (r, ri) {
    var trr = el("tr");
    var td0 = el("td", { class: "sticky-l" });
    var wrap = el("div", { class: "gene" });
    wrap.appendChild(el("span", { class: "rank", text: String(ri + 1) }));
    if (r.mito) { wrap.appendChild(el("span", { class: "dot", title: "Mitochondrial gene" })); }
    wrap.appendChild(el("span", { text: r.gene }));
    if (r.mito) wrap.appendChild(el("span", { class: "chip mito", text: "MT" }));
    td0.appendChild(wrap);
    td0.title = r.id;
    trr.appendChild(td0);
    trr.appendChild(el("td", { text: r.chr || "—" }));
    trr.appendChild(el("td", { html: "<b>" + esc(geneMean(r).toLocaleString(undefined,
      { maximumFractionDigits: geneUnit === "counts" ? 0 : 2 })) + "</b>" }));
    var vals = geneValues(r);
    var rowMax = Math.max.apply(null, vals.concat([0]));
    var rowMin = Math.min.apply(null, vals.concat([rowMax]));
    vals.forEach(function (v, i) {
      var td = el("td", { text: v.toLocaleString(undefined, { maximumFractionDigits: geneUnit === "counts" ? 0 : 2 }) });
      if (rowMax > 0) {
        var t = rowMax > rowMin ? (v - rowMin) / (rowMax - rowMin) : 0.5;
        td.style.background = "color-mix(in srgb, var(--accent) " +
          (6 + t * 26).toFixed(1) + "%, transparent)";
      }
      td.setAttribute("data-sample", samples[i]);
      trr.appendChild(td);
    });
    tbody.appendChild(trr);
  });
  tbl.appendChild(tbody);
  applyHighlight();
}

/* ========================================================== GENE TABLE NOTE */
// No CSV export here on purpose. This table is the top 100 genes by mean TPM, and a
// download button invites people to treat that subset as the gene list. The complete
// matrix ships with the results, so point at it by its real filename instead.
function renderGeneNote() {
  var host = $("#geneNote");
  if (!host) return;
  var g = D.top_genes;
  if (!g || !g.rows.length) { host.style.display = "none"; return; }

  var src = D.sources || {};
  function base(p) { return p ? String(p).split("/").pop() : null; }
  var tpm = base(src.tpm), counts = base(src.counts);

  var txt = "Showing the top " + g.rows.length + " genes by mean TPM — not the full gene list. " +
    "The complete gene-level matrix for every sample is included in your results folder";
  if (tpm || counts) {
    var files = [];
    if (tpm) files.push("<code>" + esc(tpm) + "</code> (TPM)");
    if (counts) files.push("<code>" + esc(counts) + "</code> (raw counts)");
    txt += ": " + files.join(" and ");
  }
  txt += ". Use those files for any downstream analysis.";

  host.innerHTML = '<span class="mark">i</span><span>' + txt + "</span>";
}

/* ============================================================ CSV EXPORT */
function downloadCSV(name, rows) {
  var csv = rows.map(function (r) {
    return r.map(function (c) {
      var s = (c === null || c === undefined) ? "" : String(c);
      return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    }).join(",");
  }).join("\n");
  var url = URL.createObjectURL(new Blob([csv], { type: "text/csv;charset=utf-8" }));
  var a = el("a", { href: url, download: name });
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
}

/* ================================================================ CHROME */
function buildCards() {
  var host = $("#cards");
  (D.cards || []).forEach(function (c) {
    host.appendChild(el("div", { class: "card stat" }, [
      el("div", { class: "label", text: c.label }),
      el("div", { class: "value", text: c.value }),
      el("div", { class: "sub", text: c.sub || "" })
    ]));
  });
}

function buildRunInfo() {
  var host = $("#runinfo");
  (D.run_info || []).forEach(function (kv) {
    host.appendChild(el("div", {}, [
      el("dt", { text: kv[0] }), el("dd", { text: kv[1] })
    ]));
  });
}

function buildNav() {
  var host = $("#nav");
  var secs = [["overview", "Overview"], ["qc", "QC metrics"], ["coverage", "Coverage"],
              ["metrics", "Sample table"], ["expression", "Top genes"], ["structure", "PCA"]];
  // The QC report section is only in the page when a report was supplied.
  if (document.getElementById("assessment")) secs.push(["assessment", "Verdicts"]);
  secs.forEach(function (s) {
    host.appendChild(el("a", { href: "#" + s[0], text: s[1], "data-sec": s[0] }));
  });
  var obs = new IntersectionObserver(function (entries) {
    entries.forEach(function (e) {
      if (!e.isIntersecting) return;
      host.querySelectorAll("a").forEach(function (a) {
        a.classList.toggle("active", a.getAttribute("data-sec") === e.target.id);
      });
    });
  }, { rootMargin: "-90px 0px -70% 0px", threshold: 0 });
  secs.forEach(function (s) { var n = document.getElementById(s[0]); if (n) obs.observe(n); });
}

function initTheme() {
  var stored = null;
  try { stored = localStorage.getItem("gdsc-dash-theme"); } catch (e) {}
  var prefersDark = window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches;
  setTheme(stored || (prefersDark ? "dark" : "light"));
  $("#themeBtn").addEventListener("click", function () {
    setTheme(document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark");
  });
}
function setTheme(t) {
  document.documentElement.setAttribute("data-theme", t);
  $("#themeIcon").textContent = t === "dark" ? "◑" : "◐";
  $("#themeLabel").textContent = t === "dark" ? "Light" : "Dark";
  try { localStorage.setItem("gdsc-dash-theme", t); } catch (e) {}
  // Every other chart paints with var(), so only the heatmap needs a repaint.
  if (document.querySelector("#distCard svg")) distPanel($("#distCard"));
}

/* ================================================================== INIT */
function init() {
  initTheme();
  buildNav();
  buildCards();
  buildRunInfo();

  // QC bar panels
  var qcHost = $("#qcCharts");
  var charted = (D.metric_defs || []).filter(function (d) {
    return d.chart && SAMPLES.some(function (s) {
      var v = M[s] ? M[s][d.key] : null;
      return v !== null && v !== undefined && isFinite(v);
    });
  });
  if (!charted.length) qcHost.appendChild(el("div", { class: "card empty", text: "No QC metrics parsed." }));
  charted.forEach(function (d) { qcHost.appendChild(barPanel(d)); });

  document.querySelectorAll("[data-sort]").forEach(function (b) {
    b.addEventListener("click", function () {
      barSort = b.getAttribute("data-sort");
      document.querySelectorAll("[data-sort]").forEach(function (o) {
        o.setAttribute("aria-pressed", String(o === b));
      });
      barPanels.forEach(function (f) { f(); });
      applyHighlight();
    });
  });

  coveragePanel($("#covCard"));

  // Metric table
  var gf = $("#groupFilter");
  gf.appendChild(el("option", { value: "__all__", text: "All metric groups" }));
  var seen = {};
  (D.metric_defs || []).forEach(function (d) {
    if (seen[d.group]) return;
    seen[d.group] = 1;
    gf.appendChild(el("option", { value: d.group, text: d.group }));
  });
  gf.addEventListener("change", renderMetricTable);
  $("#metricSearch").addEventListener("input", renderMetricTable);
  renderMetricTable();

  $("#dlMetrics").addEventListener("click", function () {
    var defs = activeDefs();
    var rows = [["sample"].concat(defs.map(function (d) { return d.label; }))];
    SAMPLES.forEach(function (s) {
      rows.push([s].concat(defs.map(function (d) {
        var v = M[s] ? M[s][d.key] : null;
        return (v === null || v === undefined) ? "" : v;
      })));
    });
    downloadCSV("sample_qc_metrics.csv", rows);
  });

  // Gene table
  $("#geneSearch").addEventListener("input", renderGeneTable);
  document.querySelectorAll("[data-unit]").forEach(function (b) {
    b.addEventListener("click", function () {
      geneUnit = b.getAttribute("data-unit");
      document.querySelectorAll("[data-unit]").forEach(function (o) {
        o.setAttribute("aria-pressed", String(o === b));
      });
      renderGeneTable();
    });
  });
  document.querySelectorAll("[data-mito]").forEach(function (b) {
    b.addEventListener("click", function () {
      geneMito = b.getAttribute("data-mito");
      document.querySelectorAll("[data-mito]").forEach(function (o) {
        o.setAttribute("aria-pressed", String(o === b));
      });
      renderGeneTable();
    });
  });
  renderGeneTable();
  renderGeneNote();

  // Structure
  renderPCA();
  distPanel($("#distCard"));
  varPanel($("#varCard"));

  document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") { pinned = null; hovered = null; applyHighlight(); hideTip(); }
  });
}

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
else init();

})();
