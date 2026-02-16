(() => {
  "use strict";

  const sphereCanvas = document.getElementById("sphereCanvas");
  const planeCanvas = document.getElementById("planeCanvas");
  const projectionSelect = document.getElementById("projectionSelect");
  const lockNorth = document.getElementById("lockNorth");
  const editMode = document.getElementById("editMode");
  const sampleT = document.getElementById("sampleT");
  const showJacobianViz = document.getElementById("showJacobianViz");
  const jacobianGlyphScale = document.getElementById("jacobianGlyphScale");
  const planeRadius = document.getElementById("planeRadius");
  const autoFitPlane = document.getElementById("autoFitPlane");
  const planeRadiusValue = document.getElementById("planeRadiusValue");
  const alphaInput = document.getElementById("alpha");
  const duInput = document.getElementById("du");
  const dvInput = document.getElementById("dv");
  const sidebarToggle = document.getElementById("sidebarToggle");
  const analysisSidebar = document.getElementById("analysisSidebar");
  const helpOverlay = document.getElementById("helpOverlay");
  const helpClose = document.getElementById("helpClose");
  const helpTitle = document.getElementById("helpTitle");
  const helpBody = document.getElementById("helpBody");

  const distancePanel = document.getElementById("distancePanel");
  const jacobianPanel = document.getElementById("jacobianPanel");
  const covariancePanel = document.getElementById("covariancePanel");

  const EPS = 1e-9;

  const state = {
    projection: projectionSelect.value,
    lockNorth: lockNorth.checked,
    editMode: editMode.value,
    sampleT: Number(sampleT.value),
    showJacobianViz: showJacobianViz.checked,
    jacobianGlyphScale: Number(jacobianGlyphScale.value),
    planeRadius: Number(planeRadius.value),
    alpha: Number(alphaInput.value),
    du: Number(duInput.value),
    dv: Number(dvInput.value),
    camera: {
      yaw: 0.8,
      pitch: -0.25
    },
    drag: {
      active: false,
      moved: false,
      lastX: 0,
      lastY: 0
    },
    planeView: {
      centerU: 0,
      centerV: 0
    },
    planeDrag: {
      active: false,
      lastX: 0,
      lastY: 0
    },
    planeTouch: {
      mode: "",
      prevDist: 0,
      prevMidX: 0,
      prevMidY: 0
    },
    points: {
      A: normalize([0.85, 0.2, 0.49]),
      B: normalize([-0.45, 0.78, 0.43]),
      T: [0, 0, 1]
    }
  };

  const helpText = {
    projection: {
      title: "投影类型",
      body:
        "Gnomonic 保大圆为直线，但距离与面积畸变增长快；Stereographic 保角但不保面积；Orthographic 视觉直观但会压缩远离切点区域。"
    },
    lockNorth: {
      title: "切点锁定北极",
      body:
        "开启后切点 T 固定在北极 (0,0,1)；关闭后可在球面点击任意点作为切点。投影分布和奇异区域会随 T 改变。"
    },
    editMode: {
      title: "编辑模式",
      body:
        "设置点 A/B/T：在左侧球面视图点击即可落点。A-B 定义测地线与分析路径，T 定义切平面法向。"
    },
    sampleT: {
      title: "局部分析点",
      bodyHtml:
        "<p><span class='expr'>t ∈ [0,1]</span> 表示沿 <span class='expr'>A→B</span> 的大圆插值位置。</p><p>局部 <span class='expr'>J</span>、畸变、协变性都在该点计算。</p>"
    },
    showJacobianViz: {
      title: "Jacobian 几何可视化",
      body:
        "在切平面显示分析点 X' 处的 Jacobian 几何图：线性化方块、单位圆映射椭圆和主伸缩方向。"
    },
    jacobianGlyphScale: {
      title: "Jacobian 图示尺度",
      body:
        "控制 Jacobian 几何图的参考尺度。适合在不同投影和缩放下调节可读性，不改变真实数值计算。"
    },
    planeRadius: {
      title: "切平面视图半径",
      bodyHtml:
        "<p>控制显示窗口 <span class='expr'>[-R,R]²</span>。</p><p><span class='expr'>R</span> 越小越放大，越大视野越广。</p>"
    },
    alpha: {
      title: "协变性参数 α",
      bodyHtml:
        "<p>用于坐标变换 <span class='expr'>v = lon + α sin(lat)</span>。</p><p>改变 <span class='expr'>α</span> 可观察度量分量变化与 <span class='expr'>ds²</span> 不变量。</p>"
    },
    du: {
      title: "测试向量 du",
      bodyHtml:
        "<p>协变性测试向量分量：<span class='expr'>dξ = (du,dv)</span>。</p><p>与 <span class='expr'>dv</span> 一起决定比较 <span class='expr'>ds²</span> 的方向与幅度。</p>"
    },
    dv: {
      title: "测试向量 dv",
      bodyHtml:
        "<p>协变性测试向量另一分量，和 <span class='expr'>du</span> 共同构成 <span class='expr'>dξ</span>。</p><p>建议小量，避免高阶误差。</p>"
    },
    autoFitPlane: {
      title: "自动缩放",
      body:
        "根据当前投影、A/B 测地线和球面网格采样自动估计合适 R，并带安全边距，减少手动拖动半径。"
    },
    analysisSidebar: {
      title: "分析面板",
      body:
        "侧边栏集中展示距离变化、局部微分畸变、协变性与测地线残差。可用右上角按钮收起，不影响主画布交互。"
    },
    covFormula: {
      title: "广义协变性公式",
      bodyHtml:
        "<p class='expr'>ds² = g<sub>ij</sub>dx<sup>i</sup>dx<sup>j</sup></p><p class='expr'>g′<sub>ab</sub> = (∂x<sup>i</sup>/∂u<sup>a</sup>)(∂x<sup>j</sup>/∂u<sup>b</sup>)g<sub>ij</sub></p><p class='expr'>d²x<sup>k</sup>/ds² + Γ<sup>k</sup><sub>ij</sub>(dx<sup>i</sup>/ds)(dx<sup>j</sup>/ds)=0</p>"
    },
    distancePanel: {
      title: "距离与弧长变化",
      body:
        "对比球面测地距离 dS、平面欧氏距离 dP 和投影曲线弧长 L(π(AB))，观察不同投影下距离如何被放大或压缩。"
    },
    jacobianPanel: {
      title: "局部 Jacobian / 畸变",
      body:
        "使用 Jacobian J 定量描述局部伸缩；奇异值给出主方向放缩，|det(J)| 给出面积畸变，σ1/σ2 衡量角畸变。"
    },
    covariancePanel: {
      title: "协变性与测地线验证",
      bodyHtml:
        "<p>显示 <span class='expr'>g, g′, ds²</span> 误差和 Christoffel 符号，验证坐标变换前后几何量一致。</p>"
    },
    m_proj: {
      title: "投影 π",
      bodyHtml: "<p>当前球面到切平面的映射 <span class='expr'>π : S² → T<sub>T</sub>S²</span>。</p>"
    },
    m_point: {
      title: "位置参数",
      bodyHtml: "<p>当前分析/测试点在球面上的地理坐标 <span class='expr'>(lat, lon)</span>（单位：度）。</p>"
    },
    m_ds: {
      title: "dₛ(A,B)",
      bodyHtml: "<p>球面测地距离（大圆弧长）<span class='expr'>d<sub>S</sub>(A,B)</span>。</p>"
    },
    m_dp: {
      title: "dₚ(A′,B′)",
      bodyHtml: "<p>投影后平面欧氏距离 <span class='expr'>d<sub>P</sub>(A′,B′)</span>。</p>"
    },
    m_rho: {
      title: "ρ = dₚ/dₛ",
      bodyHtml: "<p>点间距离放缩比 <span class='expr'>ρ = d<sub>P</sub>/d<sub>S</sub></span>。</p>"
    },
    m_lpi: {
      title: "L(π(AB))",
      bodyHtml: "<p>曲线 <span class='expr'>AB</span> 投影后的弧长 <span class='expr'>L(π(AB))</span>。</p>"
    },
    m_lambda: {
      title: "λ = L(π(AB))/dₛ",
      bodyHtml: "<p>弧长放缩比 <span class='expr'>λ = L(π(AB))/d<sub>S</sub></span>。</p>"
    },
    m_ca: {
      title: "cA",
      bodyHtml: "<p><span class='expr'>c<sub>A</sub> = A·T</span>，衡量点 A 相对切平面的可见性/奇异性距离。</p>"
    },
    m_cb: {
      title: "cB",
      bodyHtml: "<p><span class='expr'>c<sub>B</sub> = B·T</span>，越接近奇异条件，投影越不稳定。</p>"
    },
    m_j: {
      title: "Jacobian 分量",
      bodyHtml: "<p><span class='expr'>J = ∂(u,v)/∂(ξ¹,ξ²)</span>，表示局部线性近似。</p>"
    },
    m_jtj: {
      title: "JᵀJ",
      bodyHtml: "<p><span class='expr'>JᵀJ</span> 是拉回度量，直接反映长度和角度畸变。</p>"
    },
    m_detj: {
      title: "det(J)",
      bodyHtml: "<p><span class='expr'>det(J)</span> 的绝对值对应面积缩放。</p>"
    },
    m_sigma: {
      title: "σ₁, σ₂",
      bodyHtml: "<p>主伸缩（奇异值）<span class='expr'>σ₁, σ₂</span>，给出两个主方向的放缩倍率。</p>"
    },
    m_kappa: {
      title: "κ = σ₁/σ₂",
      bodyHtml: "<p>各向异性比 <span class='expr'>κ = σ₁/σ₂</span>，越大表示角畸变越明显。</p>"
    },
    m_theta: {
      title: "θmax",
      bodyHtml: "<p>最大角偏差估计 <span class='expr'>θ<sub>max</sub></span>（度）。</p>"
    },
    m_g: {
      title: "度量分量 g, g′",
      bodyHtml: "<p><span class='expr'>g</span> 与 <span class='expr'>g′</span> 是不同坐标表示下的度量分量。</p>"
    },
    m_ds2: {
      title: "ds² 一致性",
      bodyHtml: "<p>比较 <span class='expr'>ds²<sub>old</sub></span> 与 <span class='expr'>ds²<sub>new</sub></span>，检验协变性。</p>"
    },
    m_gamma: {
      title: "Γ 符号",
      bodyHtml: "<p>Christoffel 符号 <span class='expr'>Γ</span> 描述联络，出现在测地线方程中。</p>"
    },
    m_res: {
      title: "测地线残差",
      bodyHtml: "<p>方程残差越小，说明数值轨迹越接近理论测地线。</p>"
    }
  };

  function clamp(x, a, b) {
    return Math.max(a, Math.min(b, x));
  }

  function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  function cross(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]
    ];
  }

  function add(a, b) {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  }

  function sub(a, b) {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  }

  function scale(v, k) {
    return [v[0] * k, v[1] * k, v[2] * k];
  }

  function norm(v) {
    return Math.hypot(v[0], v[1], v[2]);
  }

  function normalize(v) {
    const n = norm(v);
    if (n < EPS) {
      return [1, 0, 0];
    }
    return [v[0] / n, v[1] / n, v[2] / n];
  }

  function fmt4(x) {
    if (!Number.isFinite(x)) {
      return "—";
    }
    return x.toFixed(4);
  }

  function fmtDeg4(rad) {
    if (!Number.isFinite(rad)) {
      return "—";
    }
    return (rad * 180 / Math.PI).toFixed(4);
  }

  function yesNo(v) {
    return v ? "是" : "否";
  }

  function renderMetricRows(el, rows) {
    el.innerHTML = rows
      .map((row) => {
        const extra = row.note ? " metric-note" : "";
        const helpBtn = row.helpKey
          ? `<button class="help-btn metric-help-btn" type="button" data-help-key="${row.helpKey}" aria-label="${row.helpAria || "参数说明"}" title="${row.helpTitle || "参数说明"}">ⓘ</button>`
          : "";
        return `<div class="metric-row${extra}"><div class="metric-label-wrap"><span class="metric-label">${row.label}</span>${helpBtn}</div><div class="metric-value">${row.value}</div></div>`;
      })
      .join("");
  }

  function mat2Vec2Mul(m, v) {
    return [m[0][0] * v[0] + m[0][1] * v[1], m[1][0] * v[0] + m[1][1] * v[1]];
  }

  function norm2(v) {
    return Math.hypot(v[0], v[1]);
  }

  function normalize2(v) {
    const n = norm2(v);
    if (n < EPS) {
      return [1, 0];
    }
    return [v[0] / n, v[1] / n];
  }

  function eigenSymmetric2x2(c11, c12, c22) {
    const tr = c11 + c22;
    const det = c11 * c22 - c12 * c12;
    const disc = Math.max(0, tr * tr - 4 * det);
    const l1 = 0.5 * (tr + Math.sqrt(disc));
    const l2 = 0.5 * (tr - Math.sqrt(disc));

    function eigenVectorFor(lambda) {
      let v = [c12, lambda - c11];
      if (Math.abs(v[0]) + Math.abs(v[1]) < 1e-10) {
        v = [lambda - c22, c12];
      }
      return normalize2(v);
    }

    return {
      lambda1: l1,
      lambda2: l2,
      v1: eigenVectorFor(l1),
      v2: eigenVectorFor(l2)
    };
  }

  function quantileSorted(sorted, q) {
    if (!sorted.length) {
      return NaN;
    }
    const t = clamp(q, 0, 1) * (sorted.length - 1);
    const i = Math.floor(t);
    const f = t - i;
    if (i + 1 >= sorted.length) {
      return sorted[i];
    }
    return sorted[i] * (1 - f) + sorted[i + 1] * f;
  }

  function setSidebarCollapsed(collapsed) {
    document.body.classList.toggle("sidebar-collapsed", collapsed);
    sidebarToggle.textContent = collapsed ? "▸" : "◂";
    const tip = collapsed ? "展开分析栏" : "收起分析栏";
    sidebarToggle.setAttribute("title", tip);
    sidebarToggle.setAttribute("aria-label", tip);
    sidebarToggle.setAttribute("aria-expanded", String(!collapsed));
    analysisSidebar.setAttribute("aria-hidden", String(collapsed));
  }

  function openHelpDialog(key) {
    const item = helpText[key] || { title: "说明", bodyHtml: "暂无说明。" };
    helpTitle.textContent = item.title;
    helpBody.innerHTML = item.bodyHtml || item.body || "暂无说明。";
    helpOverlay.hidden = false;
  }

  function closeHelpDialog() {
    helpOverlay.hidden = true;
  }

  function syncPlaneRadiusUI() {
    planeRadius.value = state.planeRadius.toFixed(1);
    planeRadiusValue.textContent = `R = ${state.planeRadius.toFixed(2)}`;
  }

  function autoFitPlaneRadius() {
    const t = state.points.T;
    const basis = tangentBasisFromNormal(t);
    const key = [];
    const context = [];

    function collectRadius(p, bucket) {
      const pr = projectPointToPlane(p, state.projection, t, basis);
      if (!pr.valid) {
        return;
      }
      const m = Math.max(Math.abs(pr.u), Math.abs(pr.v));
      if (Number.isFinite(m) && m > 1e-6 && m < 1e4) {
        bucket.push(m);
      }
    }

    collectRadius(state.points.A, key);
    collectRadius(state.points.B, key);

    const arc = sampleGeodesic(state.points.A, state.points.B, 360);
    for (let i = 0; i < arc.length; i += 1) {
      collectRadius(arc[i], key);
    }

    const latN = 32;
    const lonN = 64;
    for (let i = 0; i <= latN; i += 1) {
      const lat = -0.5 * Math.PI + (i * Math.PI) / latN;
      for (let j = 0; j < lonN; j += 1) {
        const lon = -Math.PI + (2 * Math.PI * j) / lonN;
        collectRadius(latLonToPoint(lat, lon), context);
      }
    }

    key.sort((a, b) => a - b);
    context.sort((a, b) => a - b);

    const keyR = quantileSorted(key, 0.92);
    const contextR = quantileSorted(context, 0.86);
    const candidate = Math.max(
      Number.isFinite(keyR) ? keyR : 0,
      Number.isFinite(contextR) ? contextR : 0,
      1.2
    );

    const minR = Number(planeRadius.min) || 0.2;
    const maxR = Number(planeRadius.max) || 30;
    state.planeRadius = clamp(candidate * 1.2, minR, maxR);
  }

  function rotateX(v, a) {
    const c = Math.cos(a);
    const s = Math.sin(a);
    return [v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]];
  }

  function rotateY(v, a) {
    const c = Math.cos(a);
    const s = Math.sin(a);
    return [c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]];
  }

  function worldToView(p) {
    return rotateY(rotateX(p, state.camera.pitch), state.camera.yaw);
  }

  function viewToWorld(p) {
    return rotateX(rotateY(p, -state.camera.yaw), -state.camera.pitch);
  }

  function tangentBasisFromNormal(n) {
    const ref = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const e1 = normalize(cross(ref, n));
    const e2 = normalize(cross(n, e1));
    return { e1, e2 };
  }

  function localBasisOnSphere(p) {
    const ref = Math.abs(p[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    const t1 = normalize(cross(ref, p));
    const t2 = normalize(cross(p, t1));
    return { t1, t2 };
  }

  function latLonToPoint(lat, lon) {
    const cl = Math.cos(lat);
    return [cl * Math.cos(lon), cl * Math.sin(lon), Math.sin(lat)];
  }

  function pointToLatLon(p) {
    const lat = Math.asin(clamp(p[2], -1, 1));
    const lon = Math.atan2(p[1], p[0]);
    return { lat, lon };
  }

  function projectPointToPlane(x, projection, t, basis) {
    const c = dot(x, t);
    let y = null;

    if (projection === "orthographic") {
      y = add(x, scale(t, 1 - c));
    } else if (projection === "gnomonic") {
      if (c <= 1e-6) {
        return { valid: false, reason: "c<=0", c };
      }
      y = scale(x, 1 / c);
    } else if (projection === "stereographic") {
      const d = 1 + c;
      if (d <= 1e-6) {
        return { valid: false, reason: "near antipode", c };
      }
      const part = scale(add(x, t), 2 / d);
      y = sub(part, t);
    }

    const yp = sub(y, t);
    const u = dot(yp, basis.e1);
    const v = dot(yp, basis.e2);
    return { valid: Number.isFinite(u) && Number.isFinite(v), u, v, c, y };
  }

  function differentialMapAtPoint(x, projection, t, basis, localBasis) {
    const c = dot(x, t);

    if (projection === "gnomonic" && c <= 1e-6) {
      return null;
    }
    if (projection === "stereographic" && 1 + c <= 1e-6) {
      return null;
    }

    function apply(v) {
      const dv = dot(v, t);
      if (projection === "orthographic") {
        return sub(v, scale(t, dv));
      }
      if (projection === "gnomonic") {
        return add(scale(v, 1 / c), scale(x, -dv / (c * c)));
      }
      const denom = 1 + c;
      return add(
        scale(v, 2 / denom),
        scale(add(x, t), (-2 * dv) / (denom * denom))
      );
    }

    const w1 = apply(localBasis.t1);
    const w2 = apply(localBasis.t2);

    const m11 = dot(w1, basis.e1);
    const m12 = dot(w2, basis.e1);
    const m21 = dot(w1, basis.e2);
    const m22 = dot(w2, basis.e2);

    const M = [
      [m11, m12],
      [m21, m22]
    ];

    const c11 = m11 * m11 + m21 * m21;
    const c12 = m11 * m12 + m21 * m22;
    const c22 = m12 * m12 + m22 * m22;
    const tr = c11 + c22;
    const detC = c11 * c22 - c12 * c12;
    const disc = Math.max(0, tr * tr - 4 * detC);
    const l1 = 0.5 * (tr + Math.sqrt(disc));
    const l2 = 0.5 * (tr - Math.sqrt(disc));
    const s1 = Math.sqrt(Math.max(0, l1));
    const s2 = Math.sqrt(Math.max(0, l2));
    const detM = m11 * m22 - m12 * m21;

    return {
      M,
      C: [
        [c11, c12],
        [c12, c22]
      ],
      singularValues: [Math.max(s1, s2), Math.min(s1, s2)],
      detM
    };
  }

  function greatCircleDistance(a, b) {
    return Math.acos(clamp(dot(a, b), -1, 1));
  }

  function slerp(a, b, t) {
    const omega = greatCircleDistance(a, b);
    if (omega < 1e-7) {
      return normalize(add(scale(a, 1 - t), scale(b, t)));
    }
    const s = Math.sin(omega);
    const k0 = Math.sin((1 - t) * omega) / s;
    const k1 = Math.sin(t * omega) / s;
    return normalize(add(scale(a, k0), scale(b, k1)));
  }

  function sampleGeodesic(a, b, n) {
    const points = [];
    for (let i = 0; i < n; i += 1) {
      const t = i / (n - 1);
      points.push(slerp(a, b, t));
    }
    return points;
  }

  function unwrapAngles(arr) {
    if (!arr.length) {
      return [];
    }
    const out = [arr[0]];
    for (let i = 1; i < arr.length; i += 1) {
      let d = arr[i] - arr[i - 1];
      while (d > Math.PI) {
        d -= 2 * Math.PI;
      }
      while (d < -Math.PI) {
        d += 2 * Math.PI;
      }
      out.push(out[i - 1] + d);
    }
    return out;
  }

  function projectWorldToSphereCanvas(p, ctx, canvas) {
    const q = worldToView(p);
    const cx = canvas.width * 0.5;
    const cy = canvas.height * 0.5;
    const r = Math.min(canvas.width, canvas.height) * 0.41;
    return {
      x: cx + q[0] * r,
      y: cy - q[1] * r,
      z: q[2],
      radius: r
    };
  }

  function sphereCanvasScreenToPoint(canvas, clientX, clientY) {
    const rect = canvas.getBoundingClientRect();
    const x = ((clientX - rect.left) / rect.width) * canvas.width;
    const y = ((clientY - rect.top) / rect.height) * canvas.height;

    const cx = canvas.width * 0.5;
    const cy = canvas.height * 0.5;
    const r = Math.min(canvas.width, canvas.height) * 0.41;

    const nx = (x - cx) / r;
    const ny = -(y - cy) / r;
    const rr = nx * nx + ny * ny;
    if (rr > 1) {
      return null;
    }

    const nz = Math.sqrt(Math.max(0, 1 - rr));
    const pView = [nx, ny, nz];
    return normalize(viewToWorld(pView));
  }

  function drawSegment2D(ctx, p0, p1, color, width) {
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.beginPath();
    ctx.moveTo(p0.x, p0.y);
    ctx.lineTo(p1.x, p1.y);
    ctx.stroke();
  }

  function drawArrow2D(ctx, p0, p1, color, width) {
    const dx = p1.x - p0.x;
    const dy = p1.y - p0.y;
    const len = Math.hypot(dx, dy);
    if (len < 1e-6) {
      return;
    }
    const ux = dx / len;
    const uy = dy / len;
    const head = Math.max(7, Math.min(13, 0.22 * len));
    const wing = head * 0.45;

    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.beginPath();
    ctx.moveTo(p0.x, p0.y);
    ctx.lineTo(p1.x, p1.y);
    ctx.stroke();

    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p1.x - ux * head - uy * wing, p1.y - uy * head + ux * wing);
    ctx.lineTo(p1.x - ux * head + uy * wing, p1.y - uy * head - ux * wing);
    ctx.closePath();
    ctx.fill();
  }

  function drawSphereView() {
    const ctx = sphereCanvas.getContext("2d");
    const w = sphereCanvas.width;
    const h = sphereCanvas.height;
    ctx.clearRect(0, 0, w, h);

    const cx = w * 0.5;
    const cy = h * 0.5;
    const r = Math.min(w, h) * 0.41;

    ctx.fillStyle = "#fffdf7";
    ctx.fillRect(0, 0, w, h);

    ctx.strokeStyle = "#d9ccbb";
    ctx.lineWidth = 1.1;
    ctx.beginPath();
    ctx.arc(cx, cy, r, 0, 2 * Math.PI);
    ctx.stroke();

    const latN = 12;
    const lonN = 24;
    const segN = 80;

    function drawPolyline3D(points, frontColor, backColor, width) {
      for (let i = 0; i < points.length - 1; i += 1) {
        const p0 = projectWorldToSphereCanvas(points[i], ctx, sphereCanvas);
        const p1 = projectWorldToSphereCanvas(points[i + 1], ctx, sphereCanvas);
        const zAvg = 0.5 * (p0.z + p1.z);
        const color = zAvg >= 0 ? frontColor : backColor;
        drawSegment2D(ctx, p0, p1, color, width);
      }
    }

    for (let i = 1; i < latN; i += 1) {
      const lat = -0.5 * Math.PI + (i * Math.PI) / latN;
      const line = [];
      for (let j = 0; j <= segN; j += 1) {
        const lon = -Math.PI + (2 * Math.PI * j) / segN;
        line.push(latLonToPoint(lat, lon));
      }
      drawPolyline3D(line, "#7e756a", "#d6cab9", 0.8);
    }

    for (let i = 0; i < lonN; i += 1) {
      const lon = -Math.PI + (2 * Math.PI * i) / lonN;
      const line = [];
      for (let j = 0; j <= segN; j += 1) {
        const lat = -0.5 * Math.PI + (Math.PI * j) / segN;
        line.push(latLonToPoint(lat, lon));
      }
      drawPolyline3D(line, "#7e756a", "#d6cab9", 0.8);
    }

    const arc = sampleGeodesic(state.points.A, state.points.B, 180);
    drawPolyline3D(arc, "#c23e23", "#e5ab9b", 2);

    function drawPoint(p, color, label) {
      const sp = projectWorldToSphereCanvas(p, ctx, sphereCanvas);
      const alpha = sp.z >= 0 ? 1 : 0.45;
      ctx.globalAlpha = alpha;
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(sp.x, sp.y, 6, 0, 2 * Math.PI);
      ctx.fill();
      ctx.globalAlpha = 1;
      ctx.fillStyle = "#1d1a16";
      ctx.font = "12px IBM Plex Mono, monospace";
      ctx.fillText(label, sp.x + 8, sp.y - 8);
    }

    drawPoint(state.points.A, "#8e2d17", "A");
    drawPoint(state.points.B, "#1e6f6a", "B");
    drawPoint(slerp(state.points.A, state.points.B, state.sampleT), "#5232a7", "X");
    drawPoint(state.points.T, "#a85b00", "T");

    ctx.fillStyle = "#5d564e";
    ctx.font = "12px IBM Plex Sans, sans-serif";
    ctx.fillText("前半球实线，后半球淡色", 12, 20);
  }

  function distortionColor(area) {
    if (!Number.isFinite(area) || area <= 0) {
      return "#a9a9a9";
    }
    const t = clamp(0.5 + 0.5 * Math.tanh(Math.log(area)), 0, 1);
    const r = Math.round(46 + t * 185);
    const g = Math.round(100 + (1 - t) * 90);
    const b = Math.round(180 - t * 130);
    return `rgb(${r}, ${g}, ${b})`;
  }

  function planeCanvasGeometry() {
    const w = planeCanvas.width;
    const h = planeCanvas.height;
    return {
      w,
      h,
      cx: w * 0.5,
      cy: h * 0.5,
      scalePx: Math.min(w, h) * 0.45
    };
  }

  function planeClientToUv(clientX, clientY, centerU = state.planeView.centerU, centerV = state.planeView.centerV, radius = state.planeRadius) {
    const rect = planeCanvas.getBoundingClientRect();
    const x = ((clientX - rect.left) / rect.width) * planeCanvas.width;
    const y = ((clientY - rect.top) / rect.height) * planeCanvas.height;
    const g = planeCanvasGeometry();
    return {
      u: centerU + ((x - g.cx) / g.scalePx) * radius,
      v: centerV - ((y - g.cy) / g.scalePx) * radius
    };
  }

  function clampPlaneRadius(r) {
    const minR = Number(planeRadius.min) || 0.2;
    const maxR = Number(planeRadius.max) || 30;
    return clamp(r, minR, maxR);
  }

  function drawPlaneView() {
    const ctx = planeCanvas.getContext("2d");
    const geom = planeCanvasGeometry();
    const { w, h, cx, cy, scalePx } = geom;
    ctx.clearRect(0, 0, w, h);

    ctx.fillStyle = "#fffdf7";
    ctx.fillRect(0, 0, w, h);
    const radius = state.planeRadius;
    const centerU = state.planeView.centerU;
    const centerV = state.planeView.centerV;

    function uvToCanvas(u, v) {
      return {
        x: cx + ((u - centerU) / radius) * scalePx,
        y: cy - ((v - centerV) / radius) * scalePx
      };
    }

    const t = state.points.T;
    const basis = tangentBasisFromNormal(t);

    ctx.strokeStyle = "#e2d7c8";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, cy);
    ctx.lineTo(w, cy);
    ctx.moveTo(cx, 0);
    ctx.lineTo(cx, h);
    ctx.stroke();

    ctx.strokeStyle = "#cebfa9";
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    ctx.arc(cx, cy, scalePx, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);

    const latN = 12;
    const lonN = 24;
    const segN = 80;

    function drawProjectedPolyline(points, width) {
      for (let i = 0; i < points.length - 1; i += 1) {
        const p = points[i];
        const q = points[i + 1];
        const pp = projectPointToPlane(p, state.projection, t, basis);
        const qq = projectPointToPlane(q, state.projection, t, basis);
        if (!pp.valid || !qq.valid) {
          continue;
        }
        const maxUv = Math.max(
          Math.abs(pp.u - centerU),
          Math.abs(pp.v - centerV),
          Math.abs(qq.u - centerU),
          Math.abs(qq.v - centerV)
        );
        if (maxUv > radius * 1.5) {
          continue;
        }
        const mid = normalize(add(p, q));
        const local = differentialMapAtPoint(mid, state.projection, t, basis, localBasisOnSphere(mid));
        const area = local ? Math.abs(local.detM) : NaN;
        const color = distortionColor(area);
        const c0 = uvToCanvas(pp.u, pp.v);
        const c1 = uvToCanvas(qq.u, qq.v);
        drawSegment2D(ctx, c0, c1, color, width);
      }
    }

    for (let i = 1; i < latN; i += 1) {
      const lat = -0.5 * Math.PI + (i * Math.PI) / latN;
      const line = [];
      for (let j = 0; j <= segN; j += 1) {
        const lon = -Math.PI + (2 * Math.PI * j) / segN;
        line.push(latLonToPoint(lat, lon));
      }
      drawProjectedPolyline(line, 1.1);
    }

    for (let i = 0; i < lonN; i += 1) {
      const lon = -Math.PI + (2 * Math.PI * i) / lonN;
      const line = [];
      for (let j = 0; j <= segN; j += 1) {
        const lat = -0.5 * Math.PI + (Math.PI * j) / segN;
        line.push(latLonToPoint(lat, lon));
      }
      drawProjectedPolyline(line, 1.1);
    }

    const arc = sampleGeodesic(state.points.A, state.points.B, 300);
    ctx.lineWidth = 2.3;
    ctx.strokeStyle = "#c23e23";
    ctx.beginPath();
    let hasStarted = false;
    for (let i = 0; i < arc.length; i += 1) {
      const pr = projectPointToPlane(arc[i], state.projection, t, basis);
      if (!pr.valid || Math.max(Math.abs(pr.u - centerU), Math.abs(pr.v - centerV)) > radius * 1.5) {
        hasStarted = false;
        continue;
      }
      const cpt = uvToCanvas(pr.u, pr.v);
      if (!hasStarted) {
        ctx.moveTo(cpt.x, cpt.y);
        hasStarted = true;
      } else {
        ctx.lineTo(cpt.x, cpt.y);
      }
    }
    ctx.stroke();

    function drawProjectedPoint(p, color, label) {
      const pr = projectPointToPlane(p, state.projection, t, basis);
      if (!pr.valid) {
        return;
      }
      if (Math.max(Math.abs(pr.u - centerU), Math.abs(pr.v - centerV)) > radius * 1.5) {
        return;
      }
      const cpt = uvToCanvas(pr.u, pr.v);
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(cpt.x, cpt.y, 5.7, 0, 2 * Math.PI);
      ctx.fill();
      ctx.fillStyle = "#1d1a16";
      ctx.font = "12px IBM Plex Mono, monospace";
      ctx.fillText(label, cpt.x + 8, cpt.y - 7);
    }

    function drawJacobianGlyph() {
      if (!state.showJacobianViz) {
        return;
      }

      const x = slerp(state.points.A, state.points.B, state.sampleT);
      const px = projectPointToPlane(x, state.projection, t, basis);
      if (!px.valid) {
        return;
      }

      const local = differentialMapAtPoint(x, state.projection, t, basis, localBasisOnSphere(x));
      if (!local) {
        return;
      }

      const origin = uvToCanvas(px.u, px.v);
      const M = local.M;
      const C = local.C;
      const ev = eigenSymmetric2x2(C[0][0], C[0][1], C[1][1]);
      const sigma1 = Math.sqrt(Math.max(0, ev.lambda1));
      const sigma2 = Math.sqrt(Math.max(0, ev.lambda2));
      const baseScale = clamp(state.jacobianGlyphScale, 0.05, 0.7);

      const squareCorners = [
        [-baseScale, -baseScale],
        [baseScale, -baseScale],
        [baseScale, baseScale],
        [-baseScale, baseScale]
      ];

      ctx.fillStyle = "#6e55b533";
      ctx.strokeStyle = "#5f42aa";
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      for (let i = 0; i < squareCorners.length; i += 1) {
        const mapped = mat2Vec2Mul(M, squareCorners[i]);
        const cpt = uvToCanvas(px.u + mapped[0], px.v + mapped[1]);
        if (i === 0) {
          ctx.moveTo(cpt.x, cpt.y);
        } else {
          ctx.lineTo(cpt.x, cpt.y);
        }
      }
      ctx.closePath();
      ctx.fill();
      ctx.stroke();

      ctx.strokeStyle = "#472e97";
      ctx.lineWidth = 1.8;
      ctx.beginPath();
      const segN = 72;
      for (let i = 0; i <= segN; i += 1) {
        const a = (2 * Math.PI * i) / segN;
        const mapped = mat2Vec2Mul(M, [baseScale * Math.cos(a), baseScale * Math.sin(a)]);
        const cpt = uvToCanvas(px.u + mapped[0], px.v + mapped[1]);
        if (i === 0) {
          ctx.moveTo(cpt.x, cpt.y);
        } else {
          ctx.lineTo(cpt.x, cpt.y);
        }
      }
      ctx.stroke();

      let dir1 = mat2Vec2Mul(M, ev.v1);
      let dir2 = mat2Vec2Mul(M, ev.v2);
      if (norm2(dir1) < 1e-8) {
        dir1 = [1, 0];
      }
      if (norm2(dir2) < 1e-8) {
        dir2 = [0, 1];
      }
      dir1 = normalize2(dir1);
      dir2 = normalize2(dir2);

      const axisMax = radius * 0.35;
      const a1 = clamp(baseScale * sigma1, baseScale * 0.4, axisMax);
      const a2 = clamp(baseScale * sigma2, baseScale * 0.35, axisMax);

      const axis1End = uvToCanvas(px.u + dir1[0] * a1, px.v + dir1[1] * a1);
      const axis2End = uvToCanvas(px.u + dir2[0] * a2, px.v + dir2[1] * a2);
      const axis1Neg = uvToCanvas(px.u - dir1[0] * a1, px.v - dir1[1] * a1);
      const axis2Neg = uvToCanvas(px.u - dir2[0] * a2, px.v - dir2[1] * a2);

      drawArrow2D(ctx, origin, axis1End, "#b02822", 2.1);
      drawArrow2D(ctx, origin, axis1Neg, "#b02822", 1.6);
      drawArrow2D(ctx, origin, axis2End, "#146f66", 2.1);
      drawArrow2D(ctx, origin, axis2Neg, "#146f66", 1.6);

      ctx.fillStyle = "#5232a7";
      ctx.beginPath();
      ctx.arc(origin.x, origin.y, 4.2, 0, 2 * Math.PI);
      ctx.fill();
      ctx.fillStyle = "#1d1a16";
      ctx.font = "12px IBM Plex Mono, monospace";
      ctx.fillText("X' (Jacobian)", origin.x + 8, origin.y - 8);
    }

    drawProjectedPoint(state.points.A, "#8e2d17", "A'");
    drawProjectedPoint(state.points.B, "#1e6f6a", "B'");
    drawProjectedPoint(state.points.T, "#a85b00", "T' (0,0)");
    drawJacobianGlyph();

    ctx.fillStyle = "#5d564e";
    ctx.font = "12px IBM Plex Sans, sans-serif";
    ctx.fillText("网格颜色=面积畸变 |det(J)|", 12, 20);
    ctx.fillText(`中心: (${fmt4(centerU)}, ${fmt4(centerV)})`, 12, 38);
    ctx.fillText(`窗口半径 R: ${radius.toFixed(2)} | 拖拽平移, 滚轮/双指缩放`, 12, 56);
    if (state.showJacobianViz) {
      ctx.fillText("紫色椭圆=J(单位圆), 红/绿轴=主伸缩方向", 12, 74);
    }
  }

  function geodesicPlaneArcLength(points, projection, t, basis) {
    let total = 0;
    let prev = null;

    for (let i = 0; i < points.length; i += 1) {
      const pr = projectPointToPlane(points[i], projection, t, basis);
      if (!pr.valid) {
        prev = null;
        continue;
      }
      if (Math.max(Math.abs(pr.u), Math.abs(pr.v)) > 1e4) {
        prev = null;
        continue;
      }
      if (prev) {
        total += Math.hypot(pr.u - prev.u, pr.v - prev.v);
      }
      prev = pr;
    }
    return total;
  }

  function buildDistancePanelData() {
    const t = state.points.T;
    const basis = tangentBasisFromNormal(t);

    const A = state.points.A;
    const B = state.points.B;
    const dSphere = greatCircleDistance(A, B);

    const pA = projectPointToPlane(A, state.projection, t, basis);
    const pB = projectPointToPlane(B, state.projection, t, basis);
    const dPlane = pA.valid && pB.valid ? Math.hypot(pA.u - pB.u, pA.v - pB.v) : NaN;

    const arc = sampleGeodesic(A, B, 320);
    const lPlane = geodesicPlaneArcLength(arc, state.projection, t, basis);

    return { dSphere, dPlane, lPlane, pA, pB };
  }

  function buildJacobianPanelData() {
    const x = slerp(state.points.A, state.points.B, state.sampleT);
    const t = state.points.T;
    const planeBasis = tangentBasisFromNormal(t);
    const localBasis = localBasisOnSphere(x);
    const local = differentialMapAtPoint(x, state.projection, t, planeBasis, localBasis);

    if (!local) {
      return { x, local: null };
    }

    const s1 = local.singularValues[0];
    const s2 = Math.max(local.singularValues[1], 1e-16);
    const anisotropy = s1 / s2;
    const angleErr = Math.asin(clamp((s1 - s2) / (s1 + s2), -1, 1)) * (180 / Math.PI);

    return {
      x,
      local,
      s1,
      s2,
      anisotropy,
      area: Math.abs(local.detM),
      angleErr
    };
  }

  function christoffelAtLat(lat) {
    const c = Math.cos(lat);
    const s = Math.sin(lat);
    if (Math.abs(c) < 1e-7) {
      return {
        G_lat_lonlon: NaN,
        G_lon_latlon: NaN
      };
    }
    return {
      G_lat_lonlon: s * c,
      G_lon_latlon: -Math.tan(lat)
    };
  }

  function geodesicResidual(A, B) {
    const omega = greatCircleDistance(A, B);
    if (omega < 1e-6) {
      return { max: 0, mean: 0, count: 0 };
    }

    const n = 320;
    const ds = omega / (n - 1);
    const samples = sampleGeodesic(A, B, n);
    const lat = [];
    const lonRaw = [];
    for (let i = 0; i < n; i += 1) {
      const ll = pointToLatLon(samples[i]);
      lat.push(ll.lat);
      lonRaw.push(ll.lon);
    }
    const lon = unwrapAngles(lonRaw);

    let maxR = 0;
    let sumR = 0;
    let cnt = 0;

    for (let i = 1; i < n - 1; i += 1) {
      const latP = (lat[i + 1] - lat[i - 1]) / (2 * ds);
      const lonP = (lon[i + 1] - lon[i - 1]) / (2 * ds);
      const latPP = (lat[i + 1] - 2 * lat[i] + lat[i - 1]) / (ds * ds);
      const lonPP = (lon[i + 1] - 2 * lon[i] + lon[i - 1]) / (ds * ds);

      const ch = christoffelAtLat(lat[i]);
      if (!Number.isFinite(ch.G_lat_lonlon) || !Number.isFinite(ch.G_lon_latlon)) {
        continue;
      }

      const rLat = latPP + ch.G_lat_lonlon * lonP * lonP;
      const rLon = lonPP + 2 * ch.G_lon_latlon * latP * lonP;
      const weighted = Math.hypot(rLat, Math.cos(lat[i]) * rLon);

      maxR = Math.max(maxR, weighted);
      sumR += weighted;
      cnt += 1;
    }

    return {
      max: maxR,
      mean: cnt > 0 ? sumR / cnt : 0,
      count: cnt
    };
  }

  function buildCovariancePanelData() {
    const x = slerp(state.points.A, state.points.B, state.sampleT);
    const ll = pointToLatLon(x);
    const lat = ll.lat;
    const lon = ll.lon;
    const alpha = state.alpha;

    const cosLat = Math.cos(lat);

    const g = [
      [1, 0],
      [0, cosLat * cosLat]
    ];

    const j = [
      [1, 0],
      [-alpha * cosLat, 1]
    ];

    const jT = [
      [j[0][0], j[1][0]],
      [j[0][1], j[1][1]]
    ];

    const temp = [
      [
        g[0][0] * j[0][0] + g[0][1] * j[1][0],
        g[0][0] * j[0][1] + g[0][1] * j[1][1]
      ],
      [
        g[1][0] * j[0][0] + g[1][1] * j[1][0],
        g[1][0] * j[0][1] + g[1][1] * j[1][1]
      ]
    ];

    const gPrime = [
      [
        jT[0][0] * temp[0][0] + jT[0][1] * temp[1][0],
        jT[0][0] * temp[0][1] + jT[0][1] * temp[1][1]
      ],
      [
        jT[1][0] * temp[0][0] + jT[1][1] * temp[1][0],
        jT[1][0] * temp[0][1] + jT[1][1] * temp[1][1]
      ]
    ];

    const dXi = [state.du, state.dv];
    const dX = [
      j[0][0] * dXi[0] + j[0][1] * dXi[1],
      j[1][0] * dXi[0] + j[1][1] * dXi[1]
    ];

    const dsOld = dX[0] * (g[0][0] * dX[0] + g[0][1] * dX[1]) + dX[1] * (g[1][0] * dX[0] + g[1][1] * dX[1]);
    const dsNew =
      dXi[0] * (gPrime[0][0] * dXi[0] + gPrime[0][1] * dXi[1]) +
      dXi[1] * (gPrime[1][0] * dXi[0] + gPrime[1][1] * dXi[1]);

    const residual = geodesicResidual(state.points.A, state.points.B);
    const ch = christoffelAtLat(lat);

    return {
      lat,
      lon,
      g,
      gPrime,
      dsOld,
      dsNew,
      dsErr: Math.abs(dsOld - dsNew),
      ch,
      residual
    };
  }

  function renderPanels() {
    const dist = buildDistancePanelData();
    const jac = buildJacobianPanelData();
    const cov = buildCovariancePanelData();

    const distanceRows = [
      { label: "投影 π", value: state.projection, helpKey: "m_proj", helpAria: "投影参数说明" },
      { label: "d<sub>S</sub>(A,B)", value: fmt4(dist.dSphere), helpKey: "m_ds", helpAria: "球面测地距离说明" },
      { label: "d<sub>P</sub>(A′,B′)", value: Number.isFinite(dist.dPlane) ? fmt4(dist.dPlane) : "不可定义", helpKey: "m_dp", helpAria: "平面欧氏距离说明" },
      { label: "ρ = d<sub>P</sub>/d<sub>S</sub>", value: Number.isFinite(dist.dPlane) ? fmt4(dist.dPlane / dist.dSphere) : "—", helpKey: "m_rho", helpAria: "距离比说明" },
      { label: "L(π(AB))", value: fmt4(dist.lPlane), helpKey: "m_lpi", helpAria: "投影曲线弧长说明" },
      { label: "λ = L(π(AB))/d<sub>S</sub>", value: fmt4(dist.lPlane / dist.dSphere), helpKey: "m_lambda", helpAria: "弧长比说明" },
      { label: "A 可投影", value: `${yesNo(dist.pA.valid)} · c<sub>A</sub>=${fmt4(dist.pA.c)}`, helpKey: "m_ca", helpAria: "A点可投影说明" },
      { label: "B 可投影", value: `${yesNo(dist.pB.valid)} · c<sub>B</sub>=${fmt4(dist.pB.c)}`, helpKey: "m_cb", helpAria: "B点可投影说明" }
    ];
    renderMetricRows(distancePanel, distanceRows);

    const xLL = pointToLatLon(jac.x);
    const jacRows = [
      { label: "分析点 t", value: fmt4(state.sampleT), helpKey: "sampleT", helpAria: "分析点参数说明" },
      { label: "lat / lon (°)", value: `${fmtDeg4(xLL.lat)} / ${fmtDeg4(xLL.lon)}`, helpKey: "m_point", helpAria: "位置参数说明" },
      { label: "几何图", value: `${state.showJacobianViz ? "开" : "关"} · 尺度 ${fmt4(state.jacobianGlyphScale)}`, helpKey: "showJacobianViz", helpAria: "几何图参数说明" }
    ];
    if (!jac.local) {
      jacRows.push({
        label: "提示",
        value: "当前点 Jacobian 不可定义（接近投影奇异位置）",
        note: true
      });
    } else {
      jacRows.push(
        { label: "J<sub>11</sub>, J<sub>12</sub>", value: `${fmt4(jac.local.M[0][0])}, ${fmt4(jac.local.M[0][1])}`, helpKey: "m_j", helpAria: "Jacobian 分量说明" },
        { label: "J<sub>21</sub>, J<sub>22</sub>", value: `${fmt4(jac.local.M[1][0])}, ${fmt4(jac.local.M[1][1])}`, helpKey: "m_j", helpAria: "Jacobian 分量说明" },
        { label: "(J<sup>T</sup>J)<sub>11</sub>, (J<sup>T</sup>J)<sub>12</sub>", value: `${fmt4(jac.local.C[0][0])}, ${fmt4(jac.local.C[0][1])}`, helpKey: "m_jtj", helpAria: "J转置J说明" },
        { label: "(J<sup>T</sup>J)<sub>22</sub>", value: fmt4(jac.local.C[1][1]), helpKey: "m_jtj", helpAria: "J转置J说明" },
        { label: "det(J)", value: fmt4(jac.local.detM), helpKey: "m_detj", helpAria: "det(J) 说明" },
        { label: "σ<sub>1</sub>, σ<sub>2</sub>", value: `${fmt4(jac.s1)}, ${fmt4(jac.s2)}`, helpKey: "m_sigma", helpAria: "主伸缩说明" },
        { label: "|det(J)|", value: fmt4(jac.area), helpKey: "m_detj", helpAria: "|det(J)| 说明" },
        { label: "κ = σ<sub>1</sub>/σ<sub>2</sub>", value: fmt4(jac.anisotropy), helpKey: "m_kappa", helpAria: "各向异性说明" },
        { label: "θ<sub>max</sub> (°)", value: fmt4(jac.angleErr), helpKey: "m_theta", helpAria: "角偏差说明" }
      );
    }
    renderMetricRows(jacobianPanel, jacRows);

    const covRows = [
      { label: "测试点 lat / lon (°)", value: `${fmtDeg4(cov.lat)} / ${fmtDeg4(cov.lon)}`, helpKey: "m_point", helpAria: "测试点说明" },
      { label: "α", value: fmt4(state.alpha), helpKey: "alpha", helpAria: "alpha 参数说明" },
      { label: "g<sub>11</sub>, g<sub>22</sub>", value: `${fmt4(cov.g[0][0])}, ${fmt4(cov.g[1][1])}`, helpKey: "m_g", helpAria: "度量说明" },
      { label: "g′<sub>11</sub>, g′<sub>12</sub>", value: `${fmt4(cov.gPrime[0][0])}, ${fmt4(cov.gPrime[0][1])}`, helpKey: "m_g", helpAria: "变换后度量说明" },
      { label: "g′<sub>22</sub>", value: fmt4(cov.gPrime[1][1]), helpKey: "m_g", helpAria: "变换后度量说明" },
      { label: "ds²<sub>old</sub>, ds²<sub>new</sub>", value: `${fmt4(cov.dsOld)}, ${fmt4(cov.dsNew)}`, helpKey: "m_ds2", helpAria: "线元一致性说明" },
      { label: "|Δds²|", value: fmt4(cov.dsErr), helpKey: "m_ds2", helpAria: "线元误差说明" },
      { label: "Γ<sup>φ</sup><sub>λλ</sub>", value: fmt4(cov.ch.G_lat_lonlon), helpKey: "m_gamma", helpAria: "Christoffel 说明" },
      { label: "Γ<sup>λ</sup><sub>φλ</sub>", value: fmt4(cov.ch.G_lon_latlon), helpKey: "m_gamma", helpAria: "Christoffel 说明" },
      { label: "测地线残差 max / mean", value: `${fmt4(cov.residual.max)} / ${fmt4(cov.residual.mean)}`, helpKey: "m_res", helpAria: "残差说明" }
    ];
    renderMetricRows(covariancePanel, covRows);
  }

  function renderAll() {
    syncPlaneRadiusUI();
    drawSphereView();
    drawPlaneView();
    renderPanels();
  }

  function handleControlChanges() {
    projectionSelect.addEventListener("change", () => {
      state.projection = projectionSelect.value;
      renderAll();
    });

    lockNorth.addEventListener("change", () => {
      state.lockNorth = lockNorth.checked;
      if (state.lockNorth) {
        state.points.T = [0, 0, 1];
      }
      renderAll();
    });

    editMode.addEventListener("change", () => {
      state.editMode = editMode.value;
    });

    sampleT.addEventListener("input", () => {
      state.sampleT = Number(sampleT.value);
      renderAll();
    });

    showJacobianViz.addEventListener("change", () => {
      state.showJacobianViz = showJacobianViz.checked;
      renderAll();
    });

    jacobianGlyphScale.addEventListener("input", () => {
      state.jacobianGlyphScale = Number(jacobianGlyphScale.value);
      renderAll();
    });

    planeRadius.addEventListener("input", () => {
      state.planeRadius = Number(planeRadius.value);
      renderAll();
    });

    autoFitPlane.addEventListener("click", () => {
      autoFitPlaneRadius();
      renderAll();
    });

    alphaInput.addEventListener("input", () => {
      state.alpha = Number(alphaInput.value);
      renderAll();
    });

    duInput.addEventListener("input", () => {
      state.du = Number(duInput.value);
      renderAll();
    });

    dvInput.addEventListener("input", () => {
      state.dv = Number(dvInput.value);
      renderAll();
    });

    sidebarToggle.addEventListener("click", () => {
      const collapsed = !document.body.classList.contains("sidebar-collapsed");
      setSidebarCollapsed(collapsed);
    });
  }

  function handleHelpButtons() {
    document.addEventListener("click", (ev) => {
      const btn = ev.target.closest("[data-help-key]");
      if (!btn) {
        return;
      }
      openHelpDialog(btn.dataset.helpKey);
    });

    helpClose.addEventListener("click", closeHelpDialog);

    helpOverlay.addEventListener("click", (ev) => {
      if (ev.target === helpOverlay) {
        closeHelpDialog();
      }
    });

    window.addEventListener("keydown", (ev) => {
      if (ev.key === "Escape" && !helpOverlay.hidden) {
        closeHelpDialog();
      }
    });
  }

  function handleSphereCanvasInteraction() {
    sphereCanvas.addEventListener("mousedown", (ev) => {
      state.drag.active = true;
      state.drag.moved = false;
      state.drag.lastX = ev.clientX;
      state.drag.lastY = ev.clientY;
    });

    window.addEventListener("mousemove", (ev) => {
      if (!state.drag.active) {
        return;
      }
      const dx = ev.clientX - state.drag.lastX;
      const dy = ev.clientY - state.drag.lastY;

      if (Math.abs(dx) + Math.abs(dy) > 1) {
        state.drag.moved = true;
      }

      state.camera.yaw += dx * 0.008;
      state.camera.pitch += dy * 0.008;
      state.camera.pitch = clamp(state.camera.pitch, -1.45, 1.45);

      state.drag.lastX = ev.clientX;
      state.drag.lastY = ev.clientY;
      renderAll();
    });

    window.addEventListener("mouseup", (ev) => {
      if (!state.drag.active) {
        return;
      }

      const wasMoved = state.drag.moved;
      state.drag.active = false;

      if (wasMoved) {
        return;
      }

      const picked = sphereCanvasScreenToPoint(sphereCanvas, ev.clientX, ev.clientY);
      if (!picked) {
        return;
      }

      if (state.editMode === "A") {
        state.points.A = picked;
      } else if (state.editMode === "B") {
        state.points.B = picked;
      } else if (state.editMode === "T") {
        if (state.lockNorth) {
          state.points.T = [0, 0, 1];
        } else {
          state.points.T = picked;
        }
      }

      renderAll();
    });
  }

  function handlePlaneCanvasInteraction() {
    planeCanvas.addEventListener("mousedown", (ev) => {
      state.planeDrag.active = true;
      state.planeDrag.lastX = ev.clientX;
      state.planeDrag.lastY = ev.clientY;
      planeCanvas.classList.add("dragging");
    });

    window.addEventListener("mousemove", (ev) => {
      if (!state.planeDrag.active) {
        return;
      }
      const prev = planeClientToUv(state.planeDrag.lastX, state.planeDrag.lastY);
      const now = planeClientToUv(ev.clientX, ev.clientY);
      state.planeView.centerU -= now.u - prev.u;
      state.planeView.centerV -= now.v - prev.v;
      state.planeDrag.lastX = ev.clientX;
      state.planeDrag.lastY = ev.clientY;
      renderAll();
    });

    window.addEventListener("mouseup", () => {
      state.planeDrag.active = false;
      planeCanvas.classList.remove("dragging");
    });

    planeCanvas.addEventListener("wheel", (ev) => {
      ev.preventDefault();
      const before = planeClientToUv(ev.clientX, ev.clientY);
      const factor = Math.exp(ev.deltaY * 0.0012);
      const newR = clampPlaneRadius(state.planeRadius * factor);
      state.planeRadius = newR;
      const after = planeClientToUv(ev.clientX, ev.clientY);
      state.planeView.centerU += before.u - after.u;
      state.planeView.centerV += before.v - after.v;
      renderAll();
    }, { passive: false });

    function touchDist(t0, t1) {
      return Math.hypot(t1.clientX - t0.clientX, t1.clientY - t0.clientY);
    }

    function touchMid(t0, t1) {
      return {
        x: 0.5 * (t0.clientX + t1.clientX),
        y: 0.5 * (t0.clientY + t1.clientY)
      };
    }

    planeCanvas.addEventListener("touchstart", (ev) => {
      if (ev.touches.length === 1) {
        const t = ev.touches[0];
        state.planeTouch.mode = "pan";
        state.planeTouch.prevMidX = t.clientX;
        state.planeTouch.prevMidY = t.clientY;
        planeCanvas.classList.add("dragging");
      } else if (ev.touches.length >= 2) {
        const t0 = ev.touches[0];
        const t1 = ev.touches[1];
        const mid = touchMid(t0, t1);
        state.planeTouch.mode = "pinch";
        state.planeTouch.prevDist = touchDist(t0, t1);
        state.planeTouch.prevMidX = mid.x;
        state.planeTouch.prevMidY = mid.y;
      }
    }, { passive: true });

    planeCanvas.addEventListener("touchmove", (ev) => {
      if (ev.touches.length === 1 && state.planeTouch.mode === "pan") {
        ev.preventDefault();
        const t = ev.touches[0];
        const prev = planeClientToUv(state.planeTouch.prevMidX, state.planeTouch.prevMidY);
        const now = planeClientToUv(t.clientX, t.clientY);
        state.planeView.centerU -= now.u - prev.u;
        state.planeView.centerV -= now.v - prev.v;
        state.planeTouch.prevMidX = t.clientX;
        state.planeTouch.prevMidY = t.clientY;
        renderAll();
      } else if (ev.touches.length >= 2) {
        ev.preventDefault();
        const t0 = ev.touches[0];
        const t1 = ev.touches[1];
        const mid = touchMid(t0, t1);
        const dist = Math.max(5, touchDist(t0, t1));

        const prevMidUV = planeClientToUv(state.planeTouch.prevMidX, state.planeTouch.prevMidY);
        const nowMidUV = planeClientToUv(mid.x, mid.y);
        state.planeView.centerU -= nowMidUV.u - prevMidUV.u;
        state.planeView.centerV -= nowMidUV.v - prevMidUV.v;

        const before = planeClientToUv(mid.x, mid.y);
        const factor = state.planeTouch.prevDist / dist;
        state.planeRadius = clampPlaneRadius(state.planeRadius * factor);
        const after = planeClientToUv(mid.x, mid.y);
        state.planeView.centerU += before.u - after.u;
        state.planeView.centerV += before.v - after.v;

        state.planeTouch.prevDist = dist;
        state.planeTouch.prevMidX = mid.x;
        state.planeTouch.prevMidY = mid.y;
        state.planeTouch.mode = "pinch";
        renderAll();
      }
    }, { passive: false });

    planeCanvas.addEventListener("touchend", (ev) => {
      if (ev.touches.length === 0) {
        state.planeTouch.mode = "";
        planeCanvas.classList.remove("dragging");
      } else if (ev.touches.length === 1) {
        state.planeTouch.mode = "pan";
        state.planeTouch.prevMidX = ev.touches[0].clientX;
        state.planeTouch.prevMidY = ev.touches[0].clientY;
        planeCanvas.classList.add("dragging");
      }
    });

    planeCanvas.addEventListener("touchcancel", () => {
      state.planeTouch.mode = "";
      planeCanvas.classList.remove("dragging");
    });
  }

  function handleResize() {
    function fitCanvas(canvas) {
      const dpr = Math.max(1, window.devicePixelRatio || 1);
      const cssW = canvas.clientWidth;
      const cssH = canvas.clientHeight;
      const w = Math.round(cssW * dpr);
      const h = Math.round(cssH * dpr);
      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
      }
    }

    const onResize = () => {
      fitCanvas(sphereCanvas);
      fitCanvas(planeCanvas);
      renderAll();
    };

    window.addEventListener("resize", onResize);
    onResize();
  }

  handleControlChanges();
  handleHelpButtons();
  setSidebarCollapsed(false);
  handleSphereCanvasInteraction();
  handlePlaneCanvasInteraction();
  handleResize();
  renderAll();
})();
