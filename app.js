(() => {
  "use strict";

  const sphereCanvas = document.getElementById("sphereCanvas");
  const planeCanvas = document.getElementById("planeCanvas");
  const projectionSelect = document.getElementById("projectionSelect");
  const lockNorth = document.getElementById("lockNorth");
  const editMode = document.getElementById("editMode");
  const sampleT = document.getElementById("sampleT");
  const planeRadius = document.getElementById("planeRadius");
  const alphaInput = document.getElementById("alpha");
  const duInput = document.getElementById("du");
  const dvInput = document.getElementById("dv");

  const distancePanel = document.getElementById("distancePanel");
  const jacobianPanel = document.getElementById("jacobianPanel");
  const covariancePanel = document.getElementById("covariancePanel");

  const EPS = 1e-9;

  const state = {
    projection: projectionSelect.value,
    lockNorth: lockNorth.checked,
    editMode: editMode.value,
    sampleT: Number(sampleT.value),
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
    points: {
      A: normalize([0.85, 0.2, 0.49]),
      B: normalize([-0.45, 0.78, 0.43]),
      T: [0, 0, 1]
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

  function mat2ToString(m) {
    return `[[${m[0][0].toExponential(4)}, ${m[0][1].toExponential(4)}],\n [${m[1][0].toExponential(4)}, ${m[1][1].toExponential(4)}]]`;
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

  function drawPlaneView() {
    const ctx = planeCanvas.getContext("2d");
    const w = planeCanvas.width;
    const h = planeCanvas.height;
    ctx.clearRect(0, 0, w, h);

    ctx.fillStyle = "#fffdf7";
    ctx.fillRect(0, 0, w, h);

    const cx = w * 0.5;
    const cy = h * 0.5;
    const scalePx = Math.min(w, h) * 0.45;
    const radius = state.planeRadius;

    function uvToCanvas(u, v) {
      return {
        x: cx + (u / radius) * scalePx,
        y: cy - (v / radius) * scalePx
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
        const maxUv = Math.max(Math.abs(pp.u), Math.abs(pp.v), Math.abs(qq.u), Math.abs(qq.v));
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
      if (!pr.valid || Math.max(Math.abs(pr.u), Math.abs(pr.v)) > radius * 1.5) {
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
      if (Math.max(Math.abs(pr.u), Math.abs(pr.v)) > radius * 1.5) {
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

    drawProjectedPoint(state.points.A, "#8e2d17", "A'");
    drawProjectedPoint(state.points.B, "#1e6f6a", "B'");
    drawProjectedPoint(state.points.T, "#a85b00", "T' (0,0)");

    ctx.fillStyle = "#5d564e";
    ctx.font = "12px IBM Plex Sans, sans-serif";
    ctx.fillText("网格颜色=面积畸变 |det(J)|", 12, 20);
    ctx.fillText(`当前平面窗口: [-${radius.toFixed(1)}, ${radius.toFixed(1)}]^2`, 12, 38);
  }

  function geodesicPlaneArcLength(points, projection, t, basis, radius) {
    let total = 0;
    let prev = null;

    for (let i = 0; i < points.length; i += 1) {
      const pr = projectPointToPlane(points[i], projection, t, basis);
      if (!pr.valid) {
        prev = null;
        continue;
      }
      if (Math.max(Math.abs(pr.u), Math.abs(pr.v)) > radius * 1.5) {
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
    const lPlane = geodesicPlaneArcLength(arc, state.projection, t, basis, state.planeRadius);

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

    const distanceLines = [];
    distanceLines.push(`投影: ${state.projection}`);
    distanceLines.push(`球面测地距离 dS(A,B): ${dist.dSphere.toFixed(8)}`);
    if (Number.isFinite(dist.dPlane)) {
      distanceLines.push(`平面欧氏距离 dP(A',B'): ${dist.dPlane.toFixed(8)}`);
      distanceLines.push(`点间距离比 dP/dS: ${(dist.dPlane / dist.dSphere).toFixed(8)}`);
    } else {
      distanceLines.push("平面欧氏距离 dP(A',B'): 不可定义（投影奇异/不可见）");
    }
    distanceLines.push(`投影后曲线弧长 L(π(AB)): ${dist.lPlane.toFixed(8)}`);
    distanceLines.push(`弧长比 L(π(AB))/dS: ${(dist.lPlane / dist.dSphere).toFixed(8)}`);
    distanceLines.push(`A 点可投影: ${dist.pA.valid ? "是" : "否"} | cA=${dist.pA.c.toFixed(6)}`);
    distanceLines.push(`B 点可投影: ${dist.pB.valid ? "是" : "否"} | cB=${dist.pB.c.toFixed(6)}`);
    distancePanel.textContent = distanceLines.join("\n");

    const xLL = pointToLatLon(jac.x);
    const jacLines = [];
    jacLines.push(`分析点 t=${state.sampleT.toFixed(2)} | lat=${(xLL.lat * 180 / Math.PI).toFixed(3)}° lon=${(xLL.lon * 180 / Math.PI).toFixed(3)}°`);
    if (!jac.local) {
      jacLines.push("当前点处 Jacobian 不可定义（接近投影奇异位置）");
    } else {
      jacLines.push(`J = ${mat2ToString(jac.local.M)}`);
      jacLines.push(`g_pullback = J^T J = ${mat2ToString(jac.local.C)}`);
      jacLines.push(`det(J): ${jac.local.detM.toExponential(6)}`);
      jacLines.push(`主伸缩 σ1, σ2: ${jac.s1.toExponential(6)}, ${jac.s2.toExponential(6)}`);
      jacLines.push(`面积畸变 |det(J)|: ${jac.area.toExponential(6)}`);
      jacLines.push(`角畸变各向异性 σ1/σ2: ${jac.anisotropy.toExponential(6)}`);
      jacLines.push(`最大角偏差上界(度): ${jac.angleErr.toFixed(6)}`);
    }
    jacobianPanel.textContent = jacLines.join("\n");

    const covLines = [];
    covLines.push(`协变性测试点 lat=${(cov.lat * 180 / Math.PI).toFixed(3)}° lon=${(cov.lon * 180 / Math.PI).toFixed(3)}°`);
    covLines.push(`原坐标度量 g(lat,lon) = ${mat2ToString(cov.g)}`);
    covLines.push(`变换: u=lat, v=lon+α sin(lat), α=${state.alpha.toFixed(2)}`);
    covLines.push(`新坐标度量 g'(u,v) = ${mat2ToString(cov.gPrime)}`);
    covLines.push(`ds^2_old: ${cov.dsOld.toExponential(8)} | ds^2_new: ${cov.dsNew.toExponential(8)}`);
    covLines.push(`|ds^2_old - ds^2_new|: ${cov.dsErr.toExponential(8)}`);
    covLines.push("Christoffel (lat/lon):");
    covLines.push(`Γ^lat_{lon lon} = ${cov.ch.G_lat_lonlon.toExponential(6)}`);
    covLines.push(`Γ^lon_{lat lon} = Γ^lon_{lon lat} = ${cov.ch.G_lon_latlon.toExponential(6)}`);
    covLines.push(`大圆测地线方程残差 max: ${cov.residual.max.toExponential(6)}`);
    covLines.push(`大圆测地线方程残差 mean: ${cov.residual.mean.toExponential(6)}`);
    covariancePanel.textContent = covLines.join("\n");
  }

  function renderAll() {
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

    planeRadius.addEventListener("input", () => {
      state.planeRadius = Number(planeRadius.value);
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
  handleSphereCanvasInteraction();
  handleResize();
  renderAll();
})();
