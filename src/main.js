class NURBSCurve {
    constructor(controlPoints, knots, weights) {
        const n = controlPoints.length; 
        this.degree = Math.min(n - 1, 3); 
        this.controlPoints = controlPoints;
        this.knots = knots;
        this.weights = weights;
    }

    // calculate the point in the curve, u in [0,1]
    evaluate(u) {
        const n = this.controlPoints.length - 1;
        const d = this.degree;

        const weightedControlPoints = this.controlPoints.map((point, i) => {
            const w = this.weights[i];
            return point.map(coord => coord * w); // [x*w, y*w, z*w]
        });

        // value of basis functions 
        const N = this.basisFunctions(u);

        let sumWeight = 0;
        let curvePoint = [0, 0, 0];
        for (let i = 0; i <= n; i++) {
            sumWeight += N[i] * this.weights[i];
            curvePoint = curvePoint.map((val, j) => val + N[i] * weightedControlPoints[i][j]);
        }

        if (sumWeight === 0) {
            console.warn('Sum of weights is zero. Check your control points and weights.');
            return [0, 0, 0]; // return a default value, not null
        }

        return curvePoint.map(coord => coord / sumWeight);
    }

    // basis functions
    basisFunctions(u) {
        const n = this.controlPoints.length - 1;
        const d = this.degree;
        const N = new Array(n + 1).fill(0);

        for (let i = 0; i <= n; i++) {
            N[i] = (u >= this.knots[i] && u < this.knots[i + 1]) ? 1 : 0;
        }

        for (let k = 1; k <= d; k++) {
            for (let i = 0; i <= n - k; i++) {
                const leftDenom = this.knots[i + k] - this.knots[i];
                const rightDenom = this.knots[i + k + 1] - this.knots[i + 1];
                const leftTerm = leftDenom === 0 ? 0 : ((u - this.knots[i]) / leftDenom) * N[i];
                const rightTerm = rightDenom === 0 ? 0 : ((this.knots[i + k + 1] - u) / rightDenom) * N[i + 1];
                N[i] = leftTerm + rightTerm;
            }
        }

        return N;
    }
}

// ViewMatrix to quaternion
function matrixToQuaternion(matrix) {
    const m = matrix;
    const trace = m[0] + m[5] + m[10];
    let q = [0, 0, 0, 1];

    if (trace > 0) {
        const s = Math.sqrt(trace + 1.0) * 2; // S=4*qw
        q[3] = 0.25 * s;
        q[0] = (m[9] - m[6]) / s; // qx
        q[1] = (m[2] - m[8]) / s; // qy
        q[2] = (m[4] - m[1]) / s; // qz
    } else if ((m[0] > m[5]) && (m[0] > m[10])) {
        const s = Math.sqrt(1.0 + m[0] - m[5] - m[10]) * 2; // S=4*qx
        q[3] = (m[9] - m[6]) / s;
        q[0] = 0.25 * s;
        q[1] = (m[1] + m[4]) / s;
        q[2] = (m[2] + m[8]) / s;
    } else if (m[5] > m[10]) {
        const s = Math.sqrt(1.0 + m[5] - m[0] - m[10]) * 2; // S=4*qy
        q[3] = (m[2] - m[8]) / s;
        q[0] = (m[1] + m[4]) / s;
        q[1] = 0.25 * s;
        q[2] = (m[6] + m[9]) / s;
    } else {
        const s = Math.sqrt(1.0 + m[10] - m[0] - m[5]) * 2; // S=4*qz
        q[3] = (m[4] - m[1]) / s;
        q[0] = (m[2] + m[8]) / s;
        q[1] = (m[6] + m[9]) / s;
        q[2] = 0.25 * s;
    }

    return q;
}

// quaternion to new ViewMatrix (modified, falttened by column)
function quaternionToMatrix(q) {
    const [x, y, z, w] = q;
    const xx = x * x;
    const xy = x * y;
    const xz = x * z;
    const xw = x * w;
    const yy = y * y;
    const yz = y * z;
    const yw = y * w;
    const zz = z * z;
    const zw = z * w;

    // return [
    //     1 - 2 * (yy + zz), 2 * (xy - zw), 2 * (xz + yw), 0,
    //     2 * (xy + zw), 1 - 2 * (xx + zz), 2 * (yz - xw), 0,
    //     2 * (xz - yw), 2 * (yz + xw), 1 - 2 * (xx + yy), 0,
    //     0, 0, 0, 1
    // ];

	return [
        1 - 2 * (yy + zz), 2 * (xy + zw), 2 * (xz - yw), 0,
        2 * (xy - zw), 1 - 2 * (xx + zz), 2 * (yz + xw), 0,
        2 * (xz + yw), 2 * (yz - xw), 1 - 2 * (xx + yy), 0,
        0, 0, 0, 1
    ];
}

// SLERP interpolations
function slerp(q1, q2, t) {
    const dot = q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2] + q1[3] * q2[3];
    const theta = Math.acos(dot) * t;
    const q2t = [q2[0] - dot * q1[0], q2[1] - dot * q1[1], q2[2] - dot * q1[2], q2[3] - dot * q1[3]];
    const norm = Math.sqrt(q2t[0] * q2t[0] + q2t[1] * q2t[1] + q2t[2] * q2t[2] + q2t[3] * q2t[3]);
    if (norm < 1e-6) return q1; 
    const q2tNormalized = q2t.map(q => q / norm);
    return [
        q1[0] * Math.cos(theta) + q2tNormalized[0] * Math.sin(theta),
        q1[1] * Math.cos(theta) + q2tNormalized[1] * Math.sin(theta),
        q1[2] * Math.cos(theta) + q2tNormalized[2] * Math.sin(theta),
        q1[3] * Math.cos(theta) + q2tNormalized[3] * Math.sin(theta)
    ];
}

// Description: Main file for the project
const settings = {
	renderingMode: "Octree", // "Octree", "Gaussian"
	renderResolution: 1.0,
	maxGaussians: 2000000,
	// scalingBeta: 2,
	bgColor: '#000000',
	speed: 0.07,
	lodLevel: dataSource.baseLODLevel,
	baseLevel: dataSource.baseLODLevel,
	maxLevel: dataSource.maxLODLevel,
	en_frustum_culling: false,
	debugDepth: false,
	sortTime: 'NaN',
	fps: 0,
	shDegree: 0,
	testingMode: false,
	fov: 90,
	fx: 367, // for testing
	fy: 367, // for testing
	poseId: 0,
	depthMax: 100, // for testing
}
if (dataParameter.depthMax) {
	settings.depthMax = dataParameter.depthMax;
}
let activeKeys = [];
const baseLevelQueue = [];

// get cameras from trajectory json
async function fetchCameras() {
    const response = await fetch(dataSource.trajectory);
    const data = await response.json();

    return data.map(item => ({
        id: item.id,
        img_name: item.img_name,
        width: item.width,
        height: item.height,
        viewMatrix: item.viewMatrix,
        fy: item.fy,
        fx: item.fx,
    }));
}

let cameras;

// fetchCameras().then(fetchedCameras => {
//     cameras = fetchedCameras;
// });

let reloadLod = false;
let update_count = 0;
// let default_status = false;

const startReloadLod = () => {
	if (!reloadLod) {
		reloadLod = true;
	}
};

function updateSaveCountDisplay() {
    const saveCountDisplay = document.getElementById('saveCountDisplay');
    saveCountDisplay.textContent = `View: ${savedViewMatrices.length}`;
}

let savedViewMatrices = [];
let trajectory = [];
let ifPlay = null;

function saveCurrentViewMatrix() {
    // 保存当前的 viewMatrix
    savedViewMatrices.push([...viewMatrix]);
    
    // 更新计数显示
    updateSaveCountDisplay();
    console.log("视角已保存！", viewMatrix);
}

function removeLastViewMatrix() {
    if (savedViewMatrices.length === 0) return;

    // 删除最后一个视角
    savedViewMatrices.pop();

    // 更新计数显示
    updateSaveCountDisplay();

    // 如果存在剩余的保存视角，跳到上一个视角
    if (savedViewMatrices.length > 0) {
        viewMatrix = savedViewMatrices[savedViewMatrices.length - 1];
        console.log("跳到上一个保存的视角：", viewMatrix);
        // 重新渲染以应用视角更新
        updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians)
    } else {
        console.log("没有保存的视角了！");
    }
}

function generateKnots(numControlPoints) {
    const knots = [];
    const n = numControlPoints;
    const p = Math.min(n - 1, 3);

    // 添加前 p + 1 个 0
    for (let i = 0; i <= p; i++) {
        knots.push(0);
    }

    // 添加中间的均匀分布的 knots
    for (let i = 1; i <= n - p - 1; i++) {
        knots.push(i / (n - p));
    }

    // 添加后 p + 1 个 1
    for (let i = 0; i <= p; i++) {
        knots.push(1);
    }

    return knots;
}
function transposeFromFlatArray(flatArray) {
    if (flatArray.length !== 16) {
        throw new Error("Input array must have exactly 16 elements.");
    }

    const transposed = new Array(16);

    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            transposed[i * 4 + j] = flatArray[j * 4 + i]; // 按列映射到行
        }
    }

    return transposed;
}
// 插值函数，包含平移和旋转
function saveViewMatricesToTrajectory() {
	const numInterpolatedPoints = 200;
    // 提取保存的视角矩阵的平移和旋转部分
    const controlPointsTranslation = savedViewMatrices.map(matrix => matrix.slice(12, 15));  // 提取平移部分
    // const controlPointsRotation = savedViewMatrices.map(matrix => matrix.slice(0, 9));  // 提取旋转部分的前3x3矩阵

    // 对平移进行插值
    const translatedKnots = generateKnots(controlPointsTranslation.length);
    const nurbsCurve = new NURBSCurve(controlPointsTranslation, translatedKnots, Array(controlPointsTranslation.length).fill(1));
	const interpolatedTranslations = [];

	for (let i = 0; i <= numInterpolatedPoints; i++) {
		const u = i / numInterpolatedPoints; // 计算当前的 u 值
		const point = nurbsCurve.evaluate(u); // 获取曲线上的点
		interpolatedTranslations.push(point); // 将点添加到结果数组中
	}

    // 对旋转进行插值
    const quaternionControlPoints = savedViewMatrices.map(matrix => (matrixToQuaternion(transposeFromFlatArray(matrix))));
    const interpolatedRotations = [];

    for (let i = 0; i < numInterpolatedPoints; i++) {
        const t = i / (numInterpolatedPoints - 1); // t 在 [0, 1] 之间
        const index = Math.floor(t * (quaternionControlPoints.length - 1));
        const nextIndex = Math.min(index + 1, quaternionControlPoints.length - 1);
        const q1 = quaternionControlPoints[index];
        const q2 = quaternionControlPoints[nextIndex];

        const slerpedQuaternion = slerp(q1, q2, t * (quaternionControlPoints.length - 1) - index);
        interpolatedRotations.push(slerpedQuaternion);
    }

    // 创建插值后的矩阵
    const interpolatedMatrices = [];
    for (let i = 0; i < numInterpolatedPoints; i++) {
        const translation = interpolatedTranslations[i];
        const rotation = interpolatedRotations[i];
        const newMatrix = Array(16).fill(0);
        newMatrix[15] = 1; // 设定齐次坐标

        // 设置旋转矩阵的部分
        const rotationMatrix = quaternionToMatrix(rotation);
        for (let j = 0; j < 16; j++) {
            newMatrix[j] = rotationMatrix[j];
        }

        // 设置平移部分
        newMatrix[12] = translation[0];
        newMatrix[13] = translation[1];
        newMatrix[14] = translation[2];

        interpolatedMatrices.push(newMatrix);
    }

    // 替换插值后的第一个和最后一个矩阵
    if (savedViewMatrices.length > 0) {
        interpolatedMatrices[0] = savedViewMatrices[0]; // 用第一个保存的矩阵替换插值后的第一个
        interpolatedMatrices[interpolatedMatrices.length - 1] = savedViewMatrices[savedViewMatrices.length - 1]; // 用最后一个保存的矩阵替换插值后的最后一个
    }

    trajectory = interpolatedMatrices;
	if (trajectory.length !== 0) {
		exportViewMatricesToJson();
		console.log("已保存相机轨迹！");
	}
}

function ClearTrajectory() {
	trajectory = [];
	savedViewMatrices = [];
	updateSaveCountDisplay();
}

 // 载入插值的 JSON 文件
 async function loadInterpolatedMatrices() {
	if (trajectory.length === 0) {
		const response = await fetch('interpolated_view_matrices.json');
		trajectory = await response.json();
		if (trajectory.length !== 0) {
			console.log("已导入相机轨迹！")
		};
	}
}

// 播放动画
function PlayTrajectory() {
	if (trajectory.length === 0) {
		console.log("无相机轨迹！");
		return;
	};
	if (ifPlay) return; // 如果动画已经在播放，直接返回

	var currentMatrixIndex = 0;


	// 开始播放动画，间隔时间可以根据需求调整
	ifPlay = setInterval(async () => {
		if (currentMatrixIndex >= trajectory.length) {
			// 动画播放结束，清除计时器
			clearInterval(ifPlay);
			ifPlay = null;
			console.log("播放结束！")
			return;
		}

		// 设置当前矩阵到 viewer
		viewMatrix = trajectory[currentMatrixIndex];
		await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians)
		
		// 前进到下一个矩阵
		currentMatrixIndex++;
	}, 100); // 100 毫秒间隔
}


function exportViewMatricesToJson() {
    // 将视角数据转换为 JSON 字符串
    const dataStr = JSON.stringify(trajectory, null, 2);
    const blob = new Blob([dataStr], { type: 'application/json' });
    const url = URL.createObjectURL(blob);

    // 创建一个临时的 <a> 元素，模拟文件下载
    const a = document.createElement('a');
    a.href = url;
    a.download = 'saved_view_matrices.json';
    a.click();

    // 释放 URL 对象
    URL.revokeObjectURL(url);
    console.log("视角数据已导出为 JSON 文件！");
}


function getProjectionMatrix(fx, fy, width, height) {
	const znear = 0.2;
	const zfar = 200;
	return [
		[(2 * fx) / width, 0, 0, 0],
		[0, -(2 * fy) / height, 0, 0],
		[0, 0, zfar / (zfar - znear), 1],
		[0, 0, -(zfar * znear) / (zfar - znear), 0],
	].flat();
}

function multiply4(a, b) {
	return [
		b[0] * a[0] + b[1] * a[4] + b[2] * a[8] + b[3] * a[12],
		b[0] * a[1] + b[1] * a[5] + b[2] * a[9] + b[3] * a[13],
		b[0] * a[2] + b[1] * a[6] + b[2] * a[10] + b[3] * a[14],
		b[0] * a[3] + b[1] * a[7] + b[2] * a[11] + b[3] * a[15],
		b[4] * a[0] + b[5] * a[4] + b[6] * a[8] + b[7] * a[12],
		b[4] * a[1] + b[5] * a[5] + b[6] * a[9] + b[7] * a[13],
		b[4] * a[2] + b[5] * a[6] + b[6] * a[10] + b[7] * a[14],
		b[4] * a[3] + b[5] * a[7] + b[6] * a[11] + b[7] * a[15],
		b[8] * a[0] + b[9] * a[4] + b[10] * a[8] + b[11] * a[12],
		b[8] * a[1] + b[9] * a[5] + b[10] * a[9] + b[11] * a[13],
		b[8] * a[2] + b[9] * a[6] + b[10] * a[10] + b[11] * a[14],
		b[8] * a[3] + b[9] * a[7] + b[10] * a[11] + b[11] * a[15],
		b[12] * a[0] + b[13] * a[4] + b[14] * a[8] + b[15] * a[12],
		b[12] * a[1] + b[13] * a[5] + b[14] * a[9] + b[15] * a[13],
		b[12] * a[2] + b[13] * a[6] + b[14] * a[10] + b[15] * a[14],
		b[12] * a[3] + b[13] * a[7] + b[14] * a[11] + b[15] * a[15],
	];
}

function invert4(a) {
	let b00 = a[0] * a[5] - a[1] * a[4];
	let b01 = a[0] * a[6] - a[2] * a[4];
	let b02 = a[0] * a[7] - a[3] * a[4];
	let b03 = a[1] * a[6] - a[2] * a[5];
	let b04 = a[1] * a[7] - a[3] * a[5];
	let b05 = a[2] * a[7] - a[3] * a[6];
	let b06 = a[8] * a[13] - a[9] * a[12];
	let b07 = a[8] * a[14] - a[10] * a[12];
	let b08 = a[8] * a[15] - a[11] * a[12];
	let b09 = a[9] * a[14] - a[10] * a[13];
	let b10 = a[9] * a[15] - a[11] * a[13];
	let b11 = a[10] * a[15] - a[11] * a[14];
	let det =
		b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
	if (!det) return null;
	return [
		(a[5] * b11 - a[6] * b10 + a[7] * b09) / det,
		(a[2] * b10 - a[1] * b11 - a[3] * b09) / det,
		(a[13] * b05 - a[14] * b04 + a[15] * b03) / det,
		(a[10] * b04 - a[9] * b05 - a[11] * b03) / det,
		(a[6] * b08 - a[4] * b11 - a[7] * b07) / det,
		(a[0] * b11 - a[2] * b08 + a[3] * b07) / det,
		(a[14] * b02 - a[12] * b05 - a[15] * b01) / det,
		(a[8] * b05 - a[10] * b02 + a[11] * b01) / det,
		(a[4] * b10 - a[5] * b08 + a[7] * b06) / det,
		(a[1] * b08 - a[0] * b10 - a[3] * b06) / det,
		(a[12] * b04 - a[13] * b02 + a[15] * b00) / det,
		(a[9] * b02 - a[8] * b04 - a[11] * b00) / det,
		(a[5] * b07 - a[4] * b09 - a[6] * b06) / det,
		(a[0] * b09 - a[1] * b07 + a[2] * b06) / det,
		(a[13] * b01 - a[12] * b03 - a[14] * b00) / det,
		(a[8] * b03 - a[9] * b01 + a[10] * b00) / det,
	];
}

function rotate4(a, rad, x, y, z) {
	let len = Math.hypot(x, y, z);
	x /= len;
	y /= len;
	z /= len;
	let s = Math.sin(rad);
	let c = Math.cos(rad);
	let t = 1 - c;
	let b00 = x * x * t + c;
	let b01 = y * x * t + z * s;
	let b02 = z * x * t - y * s;
	let b10 = x * y * t - z * s;
	let b11 = y * y * t + c;
	let b12 = z * y * t + x * s;
	let b20 = x * z * t + y * s;
	let b21 = y * z * t - x * s;
	let b22 = z * z * t + c;
	return [
		a[0] * b00 + a[4] * b01 + a[8] * b02,
		a[1] * b00 + a[5] * b01 + a[9] * b02,
		a[2] * b00 + a[6] * b01 + a[10] * b02,
		a[3] * b00 + a[7] * b01 + a[11] * b02,
		a[0] * b10 + a[4] * b11 + a[8] * b12,
		a[1] * b10 + a[5] * b11 + a[9] * b12,
		a[2] * b10 + a[6] * b11 + a[10] * b12,
		a[3] * b10 + a[7] * b11 + a[11] * b12,
		a[0] * b20 + a[4] * b21 + a[8] * b22,
		a[1] * b20 + a[5] * b21 + a[9] * b22,
		a[2] * b20 + a[6] * b21 + a[10] * b22,
		a[3] * b20 + a[7] * b21 + a[11] * b22,
		...a.slice(12, 16),
	];
}

function translate4(a, x, y, z) {
	return [
		...a.slice(0, 12),
		a[0] * x + a[4] * y + a[8] * z + a[12],
		a[1] * x + a[5] * y + a[9] * z + a[13],
		a[2] * x + a[6] * y + a[10] * z + a[14],
		a[3] * x + a[7] * y + a[11] * z + a[15],
	];
}

function fov2focal(fov, length) {
	return (length / 2) / Math.tan((fov * Math.PI) / 360);
}


function createWorker(self) {
	let buffer;
	let vertexCount = 0;
	let baseVertexCount = 0;

	let viewProj;
	// 6*4 + 4 + 4 = 8*4
	// XYZ - Position (Float32)
	// XYZ - Scale (Float32)
	// RGBA - colors (uint8)
	// IJKL - quaternion/rot (uint8)
	const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
	let lastProj = [];
	let depthIndex = new Uint32Array();
	let lastVertexCount = 0;

	var _floatView = new Float32Array(1);
	var _int32View = new Int32Array(_floatView.buffer);

	function floatToHalf(float) {
		_floatView[0] = float;
		var f = _int32View[0];

		var sign = (f >> 31) & 0x0001;
		var exp = (f >> 23) & 0x00ff;
		var frac = f & 0x007fffff;

		var newExp;
		if (exp == 0) {
			newExp = 0;
		} else if (exp < 113) {
			newExp = 0;
			frac |= 0x00800000;
			frac = frac >> (113 - exp);
			if (frac & 0x01000000) {
				newExp = 1;
				frac = 0;
			}
		} else if (exp < 142) {
			newExp = exp - 112;
		} else {
			newExp = 31;
			frac = 0;
		}

		return (sign << 15) | (newExp << 10) | (frac >> 13);
	}

	function packHalf2x16(x, y) {
		return (floatToHalf(x) | (floatToHalf(y) << 16)) >>> 0;
	}

	function generateTexture() {
		if (!buffer) return;
		const f_buffer = new Float32Array(buffer);
		const u_buffer = new Uint8Array(buffer);

		var texwidth = 1024 * 2; // Set to your desired width
		var texheight = Math.ceil((2 * vertexCount) / texwidth); // Set to your desired height
		var texdata = new Uint32Array(texwidth * texheight * 4); // 4 components per pixel (RGBA)
		var texdata_c = new Uint8Array(texdata.buffer);
		var texdata_f = new Float32Array(texdata.buffer);

		// Here we convert from a .splat file buffer into a texture
		// With a little bit more foresight perhaps this texture file
		// should have been the native format as it'd be very easy to
		// load it into webgl.
		for (let i = 0; i < vertexCount; i++) {
			// x, y, z
			texdata_f[8 * i + 0] = f_buffer[8 * i + 0];
			texdata_f[8 * i + 1] = f_buffer[8 * i + 1];
			texdata_f[8 * i + 2] = f_buffer[8 * i + 2];

			// r, g, b, a
			texdata_c[4 * (8 * i + 7) + 0] = u_buffer[32 * i + 24 + 0];
			texdata_c[4 * (8 * i + 7) + 1] = u_buffer[32 * i + 24 + 1];
			texdata_c[4 * (8 * i + 7) + 2] = u_buffer[32 * i + 24 + 2];
			texdata_c[4 * (8 * i + 7) + 3] = u_buffer[32 * i + 24 + 3];

			// quaternions
			let scale = [
				f_buffer[8 * i + 3 + 0],
				f_buffer[8 * i + 3 + 1],
				f_buffer[8 * i + 3 + 2],
			];
			let rot = [
				(u_buffer[32 * i + 28 + 0] - 128) / 128,
				(u_buffer[32 * i + 28 + 1] - 128) / 128,
				(u_buffer[32 * i + 28 + 2] - 128) / 128,
				(u_buffer[32 * i + 28 + 3] - 128) / 128,
			];

			// Compute the matrix product of S and R (M = S * R)
			const M = [
				1.0 - 2.0 * (rot[2] * rot[2] + rot[3] * rot[3]),
				2.0 * (rot[1] * rot[2] + rot[0] * rot[3]),
				2.0 * (rot[1] * rot[3] - rot[0] * rot[2]),

				2.0 * (rot[1] * rot[2] - rot[0] * rot[3]),
				1.0 - 2.0 * (rot[1] * rot[1] + rot[3] * rot[3]),
				2.0 * (rot[2] * rot[3] + rot[0] * rot[1]),

				2.0 * (rot[1] * rot[3] + rot[0] * rot[2]),
				2.0 * (rot[2] * rot[3] - rot[0] * rot[1]),
				1.0 - 2.0 * (rot[1] * rot[1] + rot[2] * rot[2]),
			].map((k, i) => k * scale[Math.floor(i / 3)]);

			const sigma = [
				M[0] * M[0] + M[3] * M[3] + M[6] * M[6],
				M[0] * M[1] + M[3] * M[4] + M[6] * M[7],
				M[0] * M[2] + M[3] * M[5] + M[6] * M[8],
				M[1] * M[1] + M[4] * M[4] + M[7] * M[7],-
				M[1] * M[2] + M[4] * M[5] + M[7] * M[8],
				M[2] * M[2] + M[5] * M[5] + M[8] * M[8],
			];

			texdata[8 * i + 4] = packHalf2x16(4 * sigma[0], 4 * sigma[1]);
			texdata[8 * i + 5] = packHalf2x16(4 * sigma[2], 4 * sigma[3]);
			texdata[8 * i + 6] = packHalf2x16(4 * sigma[4], 4 * sigma[5]);
		}

		self.postMessage({ texdata, texwidth, texheight }, [texdata.buffer]);
	}

	function runSort(viewProj) {
		if (!buffer) return;
		const f_buffer = new Float32Array(buffer);
		if (lastVertexCount == vertexCount) {
			let dot =
				lastProj[2] * viewProj[2] +
				lastProj[6] * viewProj[6] +
				lastProj[10] * viewProj[10];
			// if the new View is similar to previous one, then no need to re-sort
			if (Math.abs(dot - 1) < 0.1) {
				return;
			}
		} else {
			lastVertexCount = vertexCount;
		}

		const start = performance.now()
		let maxDepth = -Infinity;
		let minDepth = Infinity;
		let sizeList = new Int32Array(vertexCount);
		for (let i = 0; i < vertexCount; i++) {
			let depth =
				((viewProj[2] * f_buffer[8 * i + 0] +
					viewProj[6] * f_buffer[8 * i + 1] +
					viewProj[10] * f_buffer[8 * i + 2]) *
					4096) |
				0;
			sizeList[i] = depth;
			if (depth > maxDepth) maxDepth = depth;
			if (depth < minDepth) minDepth = depth;
		}

		// This is a 16 bit single-pass counting sort
		let depthInv = (256 * 256) / (maxDepth - minDepth);
		let counts0 = new Uint32Array(256 * 256);
		for (let i = 0; i < vertexCount; i++) {
			sizeList[i] = ((sizeList[i] - minDepth) * depthInv) | 0;
			counts0[sizeList[i]]++;
		}
		let starts0 = new Uint32Array(256 * 256);
		for (let i = 1; i < 256 * 256; i++)
			starts0[i] = starts0[i - 1] + counts0[i - 1];
		depthIndex = new Uint32Array(vertexCount);
		for (let i = 0; i < vertexCount; i++)
			// depthIndex[starts0[sizeList[i]]++] = i;
			depthIndex[starts0[sizeList[i]]++] = i;

		const sortTime = `${((performance.now() - start) / 1000).toFixed(3)}s`

		lastProj = viewProj;
		// depthIndex.reverse();
		generateTexture();
		self.postMessage({ depthIndex, viewProj, vertexCount, sortTime }, [
			depthIndex.buffer,
		]);
	}

	async function runSortMultiThread(viewProj) {
		if (!buffer) return; // 如果没有数据，直接返回
		const f_buffer = new Float32Array(buffer);
	
		// 如果顶点数量没有改变，且视图矩阵变化较小，则不需要重新排序
		if (lastVertexCount === vertexCount) {
			let dot =
				lastProj[2] * viewProj[2] +
				lastProj[6] * viewProj[6] +
				lastProj[10] * viewProj[10];
			if (Math.abs(dot - 1) < 0.1) return; // 如果视图矩阵变化小于阈值，直接返回
		} else {
			lastVertexCount = vertexCount; // 更新顶点数量
		}
	
		const start = performance.now();
	
		// 分块处理数据，减少主线程负担
		const chunkSize = 100000; // 每个分块的大小
		const chunks = [];
		for (let i = 0; i < vertexCount; i += chunkSize) {
			chunks.push(f_buffer.subarray(8 * i, 8 * Math.min(i + chunkSize, vertexCount)));
		}
	
		// 创建 Web Worker 并发送任务
		const workers = [];
		const promises = [];
		let globalMinDepth = Infinity;
		let globalMaxDepth = -Infinity;
		let depthLists = [];
	
		for (let i = 0; i < chunks.length; i++) {
			const worker = new Worker(		
				URL.createObjectURL(
					new Blob([DepthWorkerSource], {
						type: "application/javascript",
					}),
				),
			);
			workers.push(worker);
	
			promises.push(
				new Promise((resolve) => {
					worker.onmessage = function (event) {
						const { minDepth, maxDepth, depthList, chunkIndex } = event.data;
						globalMinDepth = Math.min(globalMinDepth, minDepth);
						globalMaxDepth = Math.max(globalMaxDepth, maxDepth);
						depthLists[chunkIndex] = depthList;
						resolve(); // 当任务完成时 resolve
					};
				})
			);
	
			worker.postMessage({ buffer: chunks[i].buffer, viewProj, chunkIndex: i });
		}
	
		// 等待所有 worker 任务完成
		await Promise.all(promises);
	
		// 合并深度列表
		const mergedDepthList = depthLists.flat();
		const depthInv = (256 * 256) / (globalMaxDepth - globalMinDepth);
		const counts = new Uint32Array(256 * 256);
		const sizeList = new Uint32Array(vertexCount);
	
		for (let i = 0; i < vertexCount; i++) {
			sizeList[i] = ((mergedDepthList[i] - globalMinDepth) * depthInv) | 0;
			counts[sizeList[i]]++;
		}
	
		// 统计计数并生成排序索引
		const starts = new Uint32Array(256 * 256);
		for (let i = 1; i < 256 * 256; i++) {
			starts[i] = starts[i - 1] + counts[i - 1];
		}
	
		depthIndex = new Uint32Array(vertexCount);
		for (let i = 0; i < vertexCount; i++) {
			depthIndex[starts[sizeList[i]]++] = i;
		}
	
		// 停止所有 Worker
		workers.forEach((worker) => worker.terminate());
	
		const sortTime = `${((performance.now() - start) / 1000).toFixed(3)}s`;
	
		// 更新最后一次的视图矩阵
		lastProj = viewProj;
	
		// 生成纹理并发送到主线程
		generateTexture();
		self.postMessage({ depthIndex, viewProj, vertexCount, sortTime }, [
			depthIndex.buffer,
		]);
	
	}

	const throttledSort = () => {
		if (!sortRunning) {
			sortRunning = true;
			let lastView = viewProj;
			runSort(lastView);
			// runSortMultiThread(lastView);

			setTimeout(() => {
				sortRunning = false;
				if (lastView !== viewProj) {
					throttledSort();
				}
			}, 0);
		}
	};

	let sortRunning;
	self.onmessage = (e) => {
		// console.log("Now main.worker")
		if (e.data.ply) {
			// console.log("ply")
			vertexCount = 0;
			runSort(viewProj);
			// runSortMultiThread(viewProj);
			buffer = processPlyBuffer(e.data.ply);
			vertexCount = Math.floor(buffer.byteLength / rowLength);
			postMessage({ buffer: buffer });
		} else if (e.data.buffer) {
			// console.log('buffer')
			buffer = e.data.buffer;
			vertexCount = e.data.vertexCount;
			// console.log("Recive Vertex Count", vertexCount);
		} else if (e.data.vertexCount) {
			// console.log('count')
			vertexCount = e.data.vertexCount;
		} else if (e.data.view) {
			// console.log('view')
			viewProj = e.data.view;
			throttledSort();
		}
	};
}

const vertexShaderSource = `
#version 300 es
precision highp float;
precision highp int;

uniform highp usampler2D u_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;

in vec2 position;
in int index;

out vec4 vColor;
out vec2 vPosition;

void main () {
    uvec4 cen = texelFetch(u_texture, ivec2((uint(index) & 0x3ffu) << 1, uint(index) >> 10), 0);
    vec4 cam = view * vec4(uintBitsToFloat(cen.xyz), 1);
    vec4 pos2d = projection * cam;

    float clip = 1.2 * pos2d.w;
    if (pos2d.z < -clip || pos2d.x < -clip || pos2d.x > clip || pos2d.y < -clip || pos2d.y > clip) {
        gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
        return;
    }

    uvec4 cov = texelFetch(u_texture, ivec2(((uint(index) & 0x3ffu) << 1) | 1u, uint(index) >> 10), 0);
    vec2 u1 = unpackHalf2x16(cov.x), u2 = unpackHalf2x16(cov.y), u3 = unpackHalf2x16(cov.z);
    mat3 Vrk = mat3(u1.x, u1.y, u2.x, u1.y, u2.y, u3.x, u2.x, u3.x, u3.y);

    mat3 J = mat3(
        focal.x / cam.z, 0., -(focal.x * cam.x) / (cam.z * cam.z), 
        0., -focal.y / cam.z, (focal.y * cam.y) / (cam.z * cam.z), 
        0., 0., 0.
    );

    mat3 T = transpose(mat3(view)) * J;
    mat3 cov2d = transpose(T) * Vrk * T;

    float mid = (cov2d[0][0] + cov2d[1][1]) / 2.0;
    float radius = length(vec2((cov2d[0][0] - cov2d[1][1]) / 2.0, cov2d[0][1]));
    float lambda1 = mid + radius, lambda2 = mid - radius;

    if(lambda2 < 0.0) return;
    vec2 diagonalVector = normalize(vec2(cov2d[0][1], lambda1 - cov2d[0][0]));
    vec2 majorAxis = min(sqrt(2.0 * lambda1), 1024.0) * diagonalVector;
    vec2 minorAxis = min(sqrt(2.0 * lambda2), 1024.0) * vec2(diagonalVector.y, -diagonalVector.x);

    vColor = clamp(pos2d.z/pos2d.w+1.0, 0.0, 1.0) * vec4((cov.w) & 0xffu, (cov.w >> 8) & 0xffu, (cov.w >> 16) & 0xffu, (cov.w >> 24) & 0xffu) / 255.0;
    vPosition = position;

    vec2 vCenter = vec2(pos2d) / pos2d.w;
    gl_Position = vec4(
        vCenter 
        + position.x * majorAxis / viewport 
        + position.y * minorAxis / viewport, 0.0, 1.0);

}
`.trim();

const fragmentShaderSource = `
#version 300 es
precision highp float;

in vec4 vColor;
in vec2 vPosition;

out vec4 fragColor;

void main () {
    float A = -dot(vPosition, vPosition);
    if (A < -4.0) discard;
    float B = exp(A) * vColor.a;

	// fragColor = vec4(vPosition.z / vPosition.w, 0.0, 0.0, 1.0);
    fragColor = vec4(B * vColor.rgb, B);
}

`.trim();
 
// Get defalut Viewmatrix
const cams = await fetch(dataSource.metaData);
const data = await cams.json();
let defaultViewMatrix = data.defaultViewMatrix;

let viewMatrix = defaultViewMatrix;
let lastViewMatrix = viewMatrix;
let projectionMatrix;
let stopReading = false;
let worker = null;
let octreeGeometry, octreeGeometryLoader;
const gaussianSplats = {}
import { loadPointCloud } from "./potree.js";

// Init settings GUI panel
let videoPlay = false;
let start = 0; 
function initGUI(resize) {
	const gui = new lil.GUI({ title: 'LetsGo Settings' })
	// gui.add(settings, 'renderingMode', ['Gaussian', 'Octree']).name('Rendering Mode')
	gui.add(settings, 'sortTime').name('Sort Time').disable().listen()
	gui.add(settings, 'fps').name('FPS').disable().listen()
	// gui.add(settings, 'renderResolution', 0.1, 1, 0.01).name('Preview Resolution')
	//     .onChange(() =>  resize())
	


	gui.add(settings, 'lodLevel', settings.baseLevel, settings.maxLevel, 1).name('LOD Level')
		.onChange(async value => {
			try {
				stopReading = true;

				if (value > settings.maxLevel) {
					value = settings.maxLevel;
				}

				if (value < settings.baseLevel) {
					value = settings.baseLevel;
				}
				settings.lodLevel = value;
				await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians)

			} catch (error) {
				throw error
			}
		})

	gui.add(settings, 'maxGaussians', 1e4, 1e7, 10000).name('Ext Gaussians')
		.onChange(async value => {
			try {
				stopReading = true;
				settings.maxGaussians = value;
				await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians)
			} catch (error) {
				throw error
			}
		})

	let controller = gui.add(settings, 'shDegree', 0, 2, 1).name('SH Degree')
		.onChange(async value => {
			try {
				stopReading = true;
				// await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians)
			} catch (error) {
				throw error
			}
		})
	
	// Disable the controller
	controller.domElement.style.pointerEvents = 'none';
	controller.domElement.style.opacity = '0.5';

	gui.add(settings, 'fov', 60, 120, 1).name('FOV')
		.onChange(async value => {
			try {
				stopReading = true;
				let fx = fov2focal(value, innerWidth);
				let fy = fov2focal(value, innerWidth);

				if (settings.testingMode) {
					// NOTE: manually set the canvas size for computing PSNR,...
					innerWidth = 734;
					innerHeight = 734;
					fx = settings.fx;
					fy = settings.fy;
				}
				
				projectionMatrix = getProjectionMatrix(
					fx,
					fy,
					innerWidth,
					innerHeight,
				);
				resize();
				await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians)
			} catch (error) {
				throw error
			}
		})


    // 创建一个容器 div 用于放置加号、减号和计数显示
    const controlsContainer = document.createElement('div');
    controlsContainer.style.display = 'flex';
    controlsContainer.style.alignItems = 'center';
    controlsContainer.style.margin = '5px';

    // 加号按钮 - 用于保存当前视角
    const addButton = document.createElement('button');
    addButton.textContent = '+';
	addButton.style.width = '30px';
    addButton.style.marginRight = '5px';
    addButton.onclick = saveCurrentViewMatrix;

    // 计数显示
    const saveCountDisplay = document.createElement('span');
    saveCountDisplay.textContent = 'View: 0';
	saveCountDisplay.style.width = '70px'; // 设置计数显示宽度
    saveCountDisplay.style.textAlign = 'center';
    saveCountDisplay.style.marginRight = '10px';
    saveCountDisplay.id = 'saveCountDisplay';

    // 减号按钮 - 用于删除最后保存的视角
    const removeButton = document.createElement('button');
    removeButton.textContent = '-';
	removeButton.style.width = '30px';
    removeButton.onclick = removeLastViewMatrix;

    // 导出按钮
    const exportButton = document.createElement('button');
    exportButton.textContent = 'Save';
    exportButton.style.margin = '5px';
	exportButton.style.width = '40px'
    exportButton.onclick = saveViewMatricesToTrajectory;

	// 清除按钮
	const claerButton = document.createElement('button');
	claerButton.textContent = 'Clear';
	claerButton.style.margin = '5px';
	claerButton.style.width = '40px'
	claerButton.onclick = ClearTrajectory;

	// 播放按钮
    const playButton = document.createElement('button');
    playButton.textContent = '▶';
    playButton.style.margin = '5px';
	playButton.style.width = '30px'
    playButton.onclick = PlayTrajectory;

    // 将加号、计数显示和减号按钮添加到容器
	controlsContainer.appendChild(playButton);
    controlsContainer.appendChild(addButton);
    controlsContainer.appendChild(removeButton);
	controlsContainer.appendChild(saveCountDisplay);
	controlsContainer.appendChild(exportButton);
	controlsContainer.appendChild(claerButton);

    // 将容器和导出按钮添加到 GUI 中
    gui.domElement.appendChild(controlsContainer);
	
	//
	let intervalId;
	// 创建一个新的 div 元素来包含按钮
	const div = document.createElement('div');
	div.style.display = 'flex'; // 设置 div 为 flex 布局
	div.style.justifyContent = 'flex-end'; // 将内容对齐到右边

	// 创建一个包含说明文案的HTML元素
	const description = document.createElement('div');
	description.innerHTML = ' Parameters:\
	<br /> LOD Level: Controls the level-of-detail of 3DGS model. Higher values result in more details \
	<br /><br /> Ext Gaussians: Controls the extra number of Gaussians. Increasing this value adds detail and complexity to the scene, but may also increase loading times and performance overhead. Adjust based on scene complexity and device capabilities. ';
	description.style.fontSize = '12px';
	description.style.color = '#fff';
	description.style.padding = '8px';
	description.style.opacity = 0.7;
	description.style.lineHeight = '150%';

	// 将HTML元素添加到GUI界面中
	const firstChildDiv = gui.domElement.querySelector('.children');
	firstChildDiv.appendChild(div);
	firstChildDiv.appendChild(description);

}
const chunksData = [];
async function main() {
	let carousel = true;

	const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
	let splatData = new Uint8Array(rowLength * 1000000);
	// const downsample =
	//     splatData.length / rowLength > 500000 ? 1 : 1 / devicePixelRatio;
	const downsample = settings.renderResolution;
	// Worker Init settings
	worker = new Worker(
		URL.createObjectURL(
			new Blob(["(", createWorker.toString(), ")(self)"], {
				type: "application/javascript",
			}),
		),
	);

	// GL setup
	const canvas = document.getElementById("canvas");

	const gl = canvas.getContext("webgl2", {
		antialias: false,
	});

	const vertexShader = gl.createShader(gl.VERTEX_SHADER);
	gl.shaderSource(vertexShader, vertexShaderSource);
	gl.compileShader(vertexShader);
	if (!gl.getShaderParameter(vertexShader, gl.COMPILE_STATUS))
		console.error(gl.getShaderInfoLog(vertexShader));

	const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
	gl.shaderSource(fragmentShader, fragmentShaderSource);
	gl.compileShader(fragmentShader);
	if (!gl.getShaderParameter(fragmentShader, gl.COMPILE_STATUS))
		console.error(gl.getShaderInfoLog(fragmentShader));

	const program = gl.createProgram();
	gl.attachShader(program, vertexShader);
	gl.attachShader(program, fragmentShader);
	gl.linkProgram(program);
	gl.useProgram(program);

	if (!gl.getProgramParameter(program, gl.LINK_STATUS))
		console.error(gl.getProgramInfoLog(program));

	// gl.disable(gl.DEPTH_TEST); // Disable depth testing
	// gl.enable(gl.DEPTH_TEST);
	// gl.depthFunc(gl.LESS); // 确保近处的物体覆盖远处的
	// Enable blending
	gl.enable(gl.BLEND);
	gl.blendFuncSeparate(
		gl.ONE_MINUS_DST_ALPHA,
		gl.ONE,
		gl.ONE_MINUS_DST_ALPHA,
		gl.ONE,
	);
	gl.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);

	const u_projection = gl.getUniformLocation(program, "projection");
	const u_viewport = gl.getUniformLocation(program, "viewport");
	const u_focal = gl.getUniformLocation(program, "focal");
	const u_view = gl.getUniformLocation(program, "view");

	// positions
	const triangleVertices = new Float32Array([-2, -2, 2, -2, 2, 2, -2, 2]);
	const vertexBuffer = gl.createBuffer();
	gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
	gl.bufferData(gl.ARRAY_BUFFER, triangleVertices, gl.STATIC_DRAW);
	const a_position = gl.getAttribLocation(program, "position");
	gl.enableVertexAttribArray(a_position);
	gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
	gl.vertexAttribPointer(a_position, 2, gl.FLOAT, false, 0, 0);

	var texture = gl.createTexture();
	gl.bindTexture(gl.TEXTURE_2D, texture);

	var u_textureLocation = gl.getUniformLocation(program, "u_texture");
	gl.uniform1i(u_textureLocation, 0);

	const indexBuffer = gl.createBuffer();
	const a_index = gl.getAttribLocation(program, "index");
	gl.enableVertexAttribArray(a_index);
	gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
	gl.vertexAttribIPointer(a_index, 1, gl.INT, false, 0, 0);
	gl.vertexAttribDivisor(a_index, 1);

	const resize = () => {
		let fx = fov2focal(settings.fov, innerWidth);
		let fy = fov2focal(settings.fov, innerWidth);

		if (settings.testingMode) {
			// NOTE: manually set the canvas size for computing PSNR,...
			innerWidth = 734;
			innerHeight = 734;
			fx = settings.fx;
			fy = settings.fy;
		}


		gl.uniform2fv(u_focal, new Float32Array([fx, fy]));


		projectionMatrix = getProjectionMatrix(
			fx, fy,
			innerWidth,
			innerHeight,
		);

		gl.uniform2fv(u_viewport, new Float32Array([innerWidth, innerHeight]));

		gl.canvas.width = Math.round(innerWidth / downsample);
		gl.canvas.height = Math.round(innerHeight / downsample);
		gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);

		gl.uniformMatrix4fv(u_projection, false, projectionMatrix);
	};

	window.addEventListener("resize", resize);
	resize();
	let indexBufferSize = 0;
	worker.onmessage = (e) => {
		if (e.data.buffer) {
			splatData = new Uint8Array(e.data.buffer);
			const blob = new Blob([splatData.buffer], {
				type: "application/octet-stream",
			});
			const link = document.createElement("a");
			link.download = "model.splat";
			link.href = URL.createObjectURL(blob);
			document.body.appendChild(link);
			link.click();
		} else if (e.data.texdata) {
			const { texdata, texwidth, texheight } = e.data;
			// console.log(texdata)
			gl.bindTexture(gl.TEXTURE_2D, texture);
			gl.texParameteri(
				gl.TEXTURE_2D,
				gl.TEXTURE_WRAP_S,
				gl.CLAMP_TO_EDGE,
			);
			gl.texParameteri(
				gl.TEXTURE_2D,
				gl.TEXTURE_WRAP_T,
				gl.CLAMP_TO_EDGE,
			);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

			gl.texImage2D(
				gl.TEXTURE_2D,
				0,
				gl.RGBA32UI,
				texwidth,
				texheight,
				0,
				gl.RGBA_INTEGER,
				gl.UNSIGNED_INT,
				texdata,
			);
			gl.activeTexture(gl.TEXTURE0);
			gl.bindTexture(gl.TEXTURE_2D, texture);
		} else if (e.data.depthIndex) {
			// console.log('depthindex')
			const { depthIndex, viewProj } = e.data;
			gl.bindBuffer(gl.ARRAY_BUFFER, indexBuffer);
			// 使用 gl.bufferSubData 替代 gl.bufferData，提高更新效率
			if (depthIndex.byteLength === indexBufferSize) {
				gl.bufferSubData(gl.ARRAY_BUFFER, 0, depthIndex);
			} else {
				// 如果大小发生变化，重新分配缓冲区
				gl.bufferData(gl.ARRAY_BUFFER, depthIndex, gl.DYNAMIC_DRAW);
				indexBufferSize = depthIndex.byteLength; // 更新缓冲区大小
			}
			vertexCount = e.data.vertexCount;
			settings.sortTime = e.data.sortTime;
		}
	};



	window.addEventListener("keydown", (e) => {
		// if (document.activeElement != document.body) return;
		carousel = false;
		if (!activeKeys.includes(e.code)) activeKeys.push(e.code);
		if (/\d/.test(e.key)) {
			camera = cameras[parseInt(e.key)];
			viewMatrix = getViewMatrix(camera);
		}
		if (e.code == "KeyV") {
			location.hash =
				"#" +
				JSON.stringify(
					viewMatrix.map((k) => Math.round(k * 100) / 100),
				);
		} else if (e.code === "KeyP") {
			viewMatrix = defaultViewMatrix;
			updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians);
			carousel = true;
		}
	});
	window.addEventListener("keyup", (e) => {
		activeKeys = activeKeys.filter((k) => k !== e.code);
		e.preventDefault();
		// startReloadDefault();
		down = false;
		startX = 0;
		startY = 0;
	});
	window.addEventListener("blur", () => {
		activeKeys = [];
	});

	window.addEventListener(
		"wheel",
		(e) => {
			carousel = false;
			e.preventDefault();
			const lineHeight = 10;
			const scale =
				e.deltaMode == 1
					? lineHeight
					: e.deltaMode == 2
						? innerHeight
						: 1;
			let inv = invert4(viewMatrix);
			if (e.shiftKey) {
				inv = translate4(
					inv,
					(e.deltaX * scale) / innerWidth,
					(e.deltaY * scale) / innerHeight,
					0,
				);
			} else if (e.ctrlKey || e.metaKey) {
				// inv = rotate4(inv,  (e.deltaX * scale) / innerWidth,  0, 0, 1);
				// inv = translate4(inv,  0, (e.deltaY * scale) / innerHeight, 0);
				// let preY = inv[13];
				inv = translate4(
					inv,
					0,
					0,
					(-10 * (e.deltaY * scale)) / innerHeight,
				);
				// inv[13] = preY;
			} else {
				let d = 4;
				inv = translate4(inv, 0, 0, d);
				inv = rotate4(inv, -(e.deltaX * scale) / innerWidth, 0, 1, 0);
				inv = rotate4(inv, (e.deltaY * scale) / innerHeight, 1, 0, 0);
				inv = translate4(inv, 0, 0, -d);
			}

			viewMatrix = invert4(inv);
		},
		{ passive: false },
	);

	let startX, startY;
	let down = false;
	canvas.addEventListener("mousedown", (e) => {
		carousel = false;
		// e.preventDefault();
		startX = e.clientX;
		startY = e.clientY;
		down = e.ctrlKey || e.metaKey ? 2 : 1;
	});

	canvas.addEventListener("contextmenu", (e) => {
		carousel = false;
		// e.preventDefault();
		startX = e.clientX;
		startY = e.clientY;
		down = 2;
	});

	canvas.addEventListener("mousemove", (e) => {
		// e.preventDefault();
		if (down == 1) {
			let inv = invert4(viewMatrix);
			let dx = (5 * (e.clientX - startX)) / innerWidth;
			let dy = (5 * (e.clientY - startY)) / innerHeight;
			let d = 4;

			inv = translate4(inv, 0, 0, d);
			inv = rotate4(inv, dx, 0, 1, 0);
			inv = rotate4(inv, -dy, 1, 0, 0);
			inv = translate4(inv, 0, 0, -d);
			viewMatrix = invert4(inv);
			// updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians);
			startX = e.clientX;
			startY = e.clientY;
		} else if (down == 2) {
			let inv = invert4(viewMatrix);
			// inv = rotateY(inv, );
			// let preY = inv[13];
			inv = translate4(
				inv,
				(-10 * (e.clientX - startX)) / innerWidth,
				0,
				(10 * (e.clientY - startY)) / innerHeight,
			);
			// inv[13] = preY;
			viewMatrix = invert4(inv);
			// updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians);
			startX = e.clientX;
			startY = e.clientY;
		}
	});
	canvas.addEventListener("mouseup", (e) => {
		startReloadLod();
		down = false;
		startX = 0;
		startY = 0;

	});

    let altX = 0,
        altY = 0;
    canvas.addEventListener(
        "touchstart",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1) {
                carousel = false;
                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
                down = 1;
            } else if (e.touches.length === 2) {
                // console.log('beep')
                carousel = false;
                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
                down = 1;
            }
        },
        { passive: false },
    );
    canvas.addEventListener(
        "touchmove",
        (e) => {
            e.preventDefault();
            if (e.touches.length === 1 && down) {
                let inv = invert4(viewMatrix);
                let dx = (4 * (e.touches[0].clientX - startX)) / innerWidth;
                let dy = (4 * (e.touches[0].clientY - startY)) / innerHeight;

                let d = 4;
                inv = translate4(inv, 0, 0, d);
                // inv = translate4(inv,  -x, -y, -z);
                // inv = translate4(inv,  x, y, z);
                inv = rotate4(inv, dx, 0, 1, 0);
                inv = rotate4(inv, -dy, 1, 0, 0);
                inv = translate4(inv, 0, 0, -d);

                viewMatrix = invert4(inv);

                startX = e.touches[0].clientX;
                startY = e.touches[0].clientY;
            } else if (e.touches.length === 2) {
                // alert('beep')
                const dtheta =
                    Math.atan2(startY - altY, startX - altX) -
                    Math.atan2(
                        e.touches[0].clientY - e.touches[1].clientY,
                        e.touches[0].clientX - e.touches[1].clientX,
                    );
                const dscale =
                    Math.hypot(startX - altX, startY - altY) /
                    Math.hypot(
                        e.touches[0].clientX - e.touches[1].clientX,
                        e.touches[0].clientY - e.touches[1].clientY,
                    );
                const dx =
                    (e.touches[0].clientX +
                        e.touches[1].clientX -
                        (startX + altX)) /
                    2;
                const dy =
                    (e.touches[0].clientY +
                        e.touches[1].clientY -
                        (startY + altY)) /
                    2;
                let inv = invert4(viewMatrix);
                // inv = translate4(inv,  0, 0, d);
                inv = rotate4(inv, dtheta, 0, 0, 1);

                inv = translate4(inv, -dx / innerWidth, -dy / innerHeight, 0);

                // let preY = inv[13];
                inv = translate4(inv, 0, 0, 3 * (1 - dscale));
                // inv[13] = preY;

                viewMatrix = invert4(inv);

                startX = e.touches[0].clientX;
                altX = e.touches[1].clientX;
                startY = e.touches[0].clientY;
                altY = e.touches[1].clientY;
            }
        },
        { passive: false },
    );
    canvas.addEventListener(
        "touchend",
        (e) => {
            e.preventDefault();
            down = false;
            startX = 0;
            startY = 0;
        },
        { passive: false },
    );

	let jumpDelta = 0;
	let vertexCount = 0;

	let lastFrame = 0;
	let avgFps = 0;


	window.addEventListener("gamepadconnected", (e) => {
		const gp = navigator.getGamepads()[e.gamepad.index];
		console.log(
			`Gamepad connected at index ${gp.index}: ${gp.id}. It has ${gp.buttons.length} buttons and ${gp.axes.length} axes.`,
		);
	});
	window.addEventListener("gamepaddisconnected", (e) => {
		console.log("Gamepad disconnected");
	});

	let leftGamepadTrigger, rightGamepadTrigger;

	function canvasSaveAsImage(filename) {
			// output images
		let width = innerWidth;
		let height = innerHeight;
		// 1. 从颜色缓冲区中读取像素数据v
		const pixels = new Uint8Array(width * height * 4);
		gl.readPixels(0, 0, width, height, gl.RGBA, gl.UNSIGNED_BYTE, pixels);

		// 2. 创建一个新的 ImageData 对象，并以反向的顺序填充其像素数据
		const flippedPixels = new Uint8ClampedArray(width * height * 4);
		for (let y = 0; y < height; y++) {
			for (let x = 0; x < width; x++) {
				for (let c = 0; c < 4; c++) {
					flippedPixels[(y * width + x) * 4 + c] = pixels[((height - y - 1) * width + x) * 4 + c];
				}
			}
		}
		const imageData = new ImageData(flippedPixels, width, height);

		// 3. 创建一个新的 canvas 元素，并将 ImageData 对象绘制到这个 canvas 上
		const canvas = document.createElement('canvas');
		canvas.width = width;
		canvas.height = height;
		const context = canvas.getContext('2d');
		context.putImageData(imageData, 0, 0);

		// 4. 将 canvas 的内容转换为一个 data URL
		// const dataUrl = canvas.toDataURL('image/jpeg');
		const dataUrl = canvas.toDataURL('image/png');

		// 5. 创建一个新的 a 元素，设置它的 href 属性为上一步中生成的 data URL
		const link = document.createElement('a');
		link.href = dataUrl;
		link.download = filename;

		// 6. 触发 a 元素的 click 事件
		link.click();
}
	let init = false;
	const frame = (now) => {
		let inv = invert4(viewMatrix);
		let shiftKey = activeKeys.includes("Shift") || activeKeys.includes("ShiftLeft") || activeKeys.includes("ShiftRight")

		if (activeKeys.includes("KeyW")) {
			if (shiftKey) {
				inv = translate4(inv, 0, -0.03, 0);
			} else {
				inv = translate4(inv, 0, 0, 0.1);
			}
		}
		if (activeKeys.includes("KeyS")) {
			if (shiftKey) {
				inv = translate4(inv, 0, 0.03, 0);
			} else {
				inv = translate4(inv, 0, 0, -0.1);
			}
		}
		if (activeKeys.includes("KeyA"))
			inv = translate4(inv, -0.03, 0, 0);
		//
		if (activeKeys.includes("KeyD"))
			inv = translate4(inv, 0.03, 0, 0);

		if (activeKeys.includes("KeyQ")) inv = rotate4(inv, 0.01, 0, 0, 1);
		if (activeKeys.includes("KeyE")) inv = rotate4(inv, -0.01, 0, 0, 1);

		if (activeKeys.includes("Space")) {
			start = Date.now();
			videoPlay = !videoPlay;
		} 

		const gamepads = navigator.getGamepads ? navigator.getGamepads() : [];
		let isJumping = activeKeys.includes("Ctrl");
		for (let gamepad of gamepads) {
			if (!gamepad) continue;

			const axisThreshold = 0.1; // Threshold to detect when the axis is intentionally moved
			const moveSpeed = 0.06;
			const rotateSpeed = 0.02;

			// Assuming the left stick controls translation (axes 0 and 1)
			if (Math.abs(gamepad.axes[0]) > axisThreshold) {
				inv = translate4(inv, moveSpeed * gamepad.axes[0], 0, 0);
				carousel = false;
			}
			if (Math.abs(gamepad.axes[1]) > axisThreshold) {
				inv = translate4(inv, 0, 0, -moveSpeed * gamepad.axes[1]);
				carousel = false;
			}
			if (gamepad.buttons[12].pressed || gamepad.buttons[13].pressed) {
				inv = translate4(inv, 0, -moveSpeed * (gamepad.buttons[12].pressed - gamepad.buttons[13].pressed), 0);
				carousel = false;
			}

			if (gamepad.buttons[14].pressed || gamepad.buttons[15].pressed) {
				inv = translate4(inv, -moveSpeed * (gamepad.buttons[14].pressed - gamepad.buttons[15].pressed), 0, 0);
				carousel = false;
			}

			// Assuming the right stick controls rotation (axes 2 and 3)
			if (Math.abs(gamepad.axes[2]) > axisThreshold) {
				inv = rotate4(inv, rotateSpeed * gamepad.axes[2], 0, 1, 0);
				carousel = false;
			}
			if (Math.abs(gamepad.axes[3]) > axisThreshold) {
				inv = rotate4(inv, -rotateSpeed * gamepad.axes[3], 1, 0, 0);
				carousel = false;
			}

			let tiltAxis = gamepad.buttons[6].value - gamepad.buttons[7].value;
			if (Math.abs(tiltAxis) > axisThreshold) {
				inv = rotate4(inv, rotateSpeed * tiltAxis, 0, 0, 1);
				carousel = false;
			}
			if (gamepad.buttons[4].pressed && !leftGamepadTrigger) {
				camera = cameras[(cameras.indexOf(camera) + 1) % cameras.length]
				inv = invert4(getViewMatrix(camera));
				carousel = false;
			}
			if (gamepad.buttons[5].pressed && !rightGamepadTrigger) {
				camera = cameras[(cameras.indexOf(camera) + cameras.length - 1) % cameras.length]
				inv = invert4(getViewMatrix(camera));
				carousel = false;
			}
			leftGamepadTrigger = gamepad.buttons[4].pressed;
			rightGamepadTrigger = gamepad.buttons[5].pressed;
			if (gamepad.buttons[0].pressed) {
				isJumping = true;
				carousel = false;
			}
			if (gamepad.buttons[3].pressed) {
				carousel = true;
			}
		}

		viewMatrix = invert4(inv);

		if (carousel) {
			let inv = invert4(defaultViewMatrix);

			const t = Math.sin((Date.now() - start) / 5000);
			inv = translate4(inv, 2.5 * t, 0, 6 * (1 - Math.cos(t)));
			inv = rotate4(inv, -0.6 * t, 0, 1, 0);

			viewMatrix = invert4(inv);
		}

		const viewProj = multiply4(projectionMatrix, viewMatrix);
		worker.postMessage({ view: viewProj });

		const lastViewProj = multiply4(projectionMatrix, lastViewMatrix);
		const viewProjNorm = Math.sqrt(viewProj[2] * viewProj[2] + viewProj[6] * viewProj[6] + viewProj[10] * viewProj[10]);
		const lastViewProjNorm = Math.sqrt(lastViewProj[2] * lastViewProj[2] + lastViewProj[6] * lastViewProj[6] + lastViewProj[10] * lastViewProj[10]);

		let dot =
			lastViewProj[2] / lastViewProjNorm * viewProj[2] / viewProjNorm +
			lastViewProj[6] / lastViewProjNorm * viewProj[6] / viewProjNorm +
			lastViewProj[10] / lastViewProjNorm * viewProj[10] / viewProjNorm;


		// if ((reloadLod || activeKeys.includes("ControlLeft")) && settings.renderingMode == "Octree" && Math.abs(dot - 1) > 0.001) {
		if ((reloadLod || activeKeys.includes("ControlLeft")) && settings.renderingMode == "Octree" ) {
			reloadLod = false;
			// console.log(viewMatrix)
			stopReading = true;
			updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians);
			lastViewMatrix = viewMatrix;
		}

		const currentFps = 1000 / ( now - lastFrame) || 0;
		avgFps = avgFps * 0.9 + currentFps * 0.1;

		if (vertexCount > 0) {
			// document.getElementById("spinner").style.display = "none";
			gl.uniformMatrix4fv(u_view, false, viewMatrix);
			gl.clear(gl.COLOR_BUFFER_BIT);
			gl.drawArraysInstanced(gl.TRIANGLE_FAN, 0, 4, vertexCount);
		} else {
			gl.clear(gl.COLOR_BUFFER_BIT);
			// document.getElementById("spinner").style.display = "";
			start = Date.now() + 2000;
		}

		settings.fps = Math.round(avgFps)
		lastFrame = now;
		requestAnimationFrame(frame);

		if (activeKeys.includes("KeyC")) {
			// image capture
			canvasSaveAsImage('screenshot.png');
		}

		if (activeKeys.includes("KeyR") && settings.testingMode) {
			// video record
			const recordCameras = async () => {
				for (let i = 0; i < cameras.length; i++) {
					viewMatrix = cameras[i].viewMatrix;
					settings.poseId = i;
					reloadLod = true;
					await updateGaussianByView(viewMatrix, projectionMatrix, settings.lodLevel, settings.maxGaussians);
					// wait 0.5 second
					await new Promise(resolve => setTimeout(resolve, 1500));
					console.log('poseId', settings.poseId);
			
					// Use requestAnimationFrame to ensure WebGL has finished rendering
					await new Promise(resolve => requestAnimationFrame(resolve));
			
					// canvasSaveAsImage(`image_${i}.jpg`);
					canvasSaveAsImage(`${String(i).padStart(5, '0')}.png`);
				}
			};

			recordCameras();
		}
	};

	frame();
	// requestAnimationFrame(frame);

	const selectFile = (file) => {
		const fr = new FileReader();

		// change render mode to gaussian-lod
		settings.mode = "gaussian";

		if (/\.json$/i.test(file.name)) {
			fr.onload = () => {
				cameras = JSON.parse(fr.result);
				viewMatrix = getViewMatrix(cameras[0]);
				fx, fy = fov2focal(settings.fov, innerWidth, innerHeight);
				projectionMatrix = getProjectionMatrix(
					fx / downsample,
					fy / downsample,
					canvas.width,
					canvas.height,
				);
				gl.uniformMatrix4fv(u_projection, false, projectionMatrix);

				console.log("Loaded Cameras");
			};
			fr.readAsText(file);
		} else {
			stopLoading = true;
			fr.onload = () => {
				splatData = new Uint8Array(fr.result);
				console.log("Loaded", Math.floor(splatData.length / rowLength));

				if (
					splatData[0] == 112 &&
					splatData[1] == 108 &&
					splatData[2] == 121 &&
					splatData[3] == 10
				) {
					// ply file magic header means it should be handled differently
					worker.postMessage({ ply: splatData.buffer });
				} else {
					worker.postMessage({
						buffer: splatData.buffer,
						vertexCount: Math.floor(splatData.length / rowLength),
					});
				}
			};
			fr.readAsArrayBuffer(file);
		}
	};

	window.addEventListener("hashchange", (e) => {
		try {
			viewMatrix = JSON.parse(decodeURIComponent(location.hash.slice(1)));
			carousel = false;
		} catch (err) { }
	});

	const preventDefault = (e) => {
		e.preventDefault();
		e.stopPropagation();
	};
	document.addEventListener("dragenter", preventDefault);
	document.addEventListener("dragover", preventDefault);

	// Setup Setting GUI
	initGUI(resize)

	// Init Gaussian Scene
	let bytesRead = 0;
	let lastVertexCount = -1;
	let stopLoading = false;

}

main().catch((err) => {
	// document.getElementById("spinner").style.display = "none";
	document.getElementById("message").innerText = err.toString();
});


function nodeView(i, node) {
	const att = node.geometry.attributes
	const min = node.boundingBox.min
	let position = [att.position.array[i * 3] + min.x,
	att.position.array[i * 3 + 1] + min.y,
	att.position.array[i * 3 + 2] + min.z]

	let harmonic = [att.f_dc_0.array[i], att.f_dc_1.array[i], att.f_dc_2.array[i]]

	let opacity = att.opacity.array[i]
	let scale = [att.scale_0.array[i], att.scale_1.array[i], att.scale_2.array[i]]
	let rotation = [att.rot_0.array[i], att.rot_1.array[i], att.rot_2.array[i], att.rot_3.array[i]]

	let deg = 0
	let sh = [att.f_dc_0.array[i], att.f_dc_1.array[i], att.f_dc_2.array[i]]

	if (att['f_rest_0'] !== undefined) {
		for (let j = 0; j <= 8; j++) {
			harmonic.push(att[`f_rest_${j}`].array[i]);
		}
		deg = 1
	}

	if (att['f_rest_9'] !== undefined) {
		for (let k = 9; k <= 23; k++) {
			harmonic.push(att[`f_rest_${k}`].array[i]);
		}
		deg = 2
	}
	// NOTE: resign the value, from GaussianView.cpp +D 2024.4.12
	for (let j = 1; j < (deg + 1) * (deg + 1); j++) {
		sh.push(harmonic[j - 1 + 3]);
		sh.push(harmonic[j - 1 + (deg + 1) * (deg + 1) + 2]);
		sh.push(harmonic[j - 1 + 2 * (deg + 1) * (deg + 1) + 1]);
	}

	harmonic = sh

	return { position, harmonic, opacity, scale, rotation }
}

async function loadOctreeGeometry(rootNode) {

	const start = performance.now()
	const sigmoid = (m1) => 1 / (1 + Math.exp(-m1))
	// Scene bouding box
	const queue = [rootNode]
	let vertexCount = 0
	let campos = [viewMatrix[2], viewMatrix[6], viewMatrix[10]];
	let add_level = true;
	// first for loop to get the count of gaussians
	while (queue.length > 0) {
		const currentNode = queue.shift();

		if (currentNode.level === settings.baseLevel) {
			const currentCount = currentNode.numPoints;
			vertexCount += currentCount;
			baseLevelQueue.push({ node: currentNode, level: currentNode.level });
			add_level = false;
		}
		// set the children of the current node

		if (currentNode.children && add_level) {
			for (let cid = 0; cid < 8; cid++) {
				const child = currentNode.children[cid];
				if (child && "geometry" in child) {
					queue.push(child);
				}
			}
		}
	}
	// 6*4 + 4 + 4 = 8*4
	// XYZ - Position (Float32)
	// XYZ - Scale (Float32)
	// RGBA - colors (uint8)
	// IJKL - quaternion/rot (uint8)
	const rowLength = 3 * 4 + 3 * 4 + 4 + 4;
	const buffer = new ArrayBuffer(rowLength * vertexCount);

	// second for loop to get the gaussian data 
	let glabol_i = 0
	queue.push(rootNode)
	while (queue.length > 0) {
		const currentNode = queue.shift();

		// set the children of the current node
		if (currentNode.children) {
			for (let cid = 0; cid < 8; cid++) {
				const child = currentNode.children[cid];
				if (child && "geometry" in child) {
					queue.push(child);
				}
			}
		}

		if (currentNode.level != settings.baseLevel) {
			continue;
		}
		if (currentNode.level > settings.baseLevel) {
			break;
		}

		const currentCount = currentNode.numPoints
		for (let i = 0; i < currentCount; i++) {
			const positions = new Float32Array(buffer, glabol_i * rowLength, 3);
			const scales = new Float32Array(buffer, glabol_i * rowLength + 4 * 3, 3);
			const rgbas = new Uint8ClampedArray(
				buffer,
				glabol_i * rowLength + 4 * 3 + 4 * 3,
				4,
			);
			const rots = new Uint8ClampedArray(
				buffer,
				glabol_i * rowLength + 4 * 3 + 4 * 3 + 4,
				4,
			);

			let { position, harmonic, opacity, scale, rotation } = nodeView(i, currentNode)
			// Normalize quaternion
			let length2 = 0

			for (let j = 0; j < 4; j++)
				length2 += rotation[j] * rotation[j]

			const length = Math.sqrt(length2)

			rotation = rotation.map(v => (v / length) * 128 + 128)

			rots[0] = rotation[0];
			rots[1] = rotation[1];
			rots[2] = rotation[2];
			rots[3] = rotation[3];

			// Exponentiate scale
			scale = scale.map(v => Math.exp(v))

			scales[0] = scale[0];
			scales[1] = scale[1];
			scales[2] = scale[2];

			let color = computeColorFromSH(settings.shDegree, position, campos, harmonic)
			rgbas[0] = color.x * 255;
			rgbas[1] = color.y * 255;
			rgbas[2] = color.z * 255;

			// Activate alpha
			rgbas[3] = sigmoid(opacity) * 255;

			positions[0] = position[0];
			positions[1] = position[1];
			positions[2] = position[2];

			glabol_i++;

		}

	}

	const loadTime = `${((performance.now() - start) / 1000).toFixed(3)}s`
	console.log(`[Loader] load ${vertexCount} gaussians in ${loadTime}.`)

	// load finished
	document.getElementById("progress").style.display = "none";
	document.getElementById("spinner").style.display = "none";
	
	return buffer;
}
let octreeFileUrl;
async function loadScene() {
	const content = await loadPointCloud(dataSource.metaData, 'geometry', octreeFileUrl)
	octreeGeometry = content.geometry
	octreeGeometryLoader = octreeGeometry.loader
	await new Promise((resolve) => setTimeout(resolve, dataSource.initLoadTime));
	// wait for the geometry to load
	console.log('INFO: octreeGeometry loading...')
	settings.renderingMode = "Octree";

	const queue = [{ node: octreeGeometry.root, level: 0 }];
	while (queue.length > 0) {
		const { node, level } = queue.shift();
		if (level <= settings.baseLevel) {
			// console.log(node)
			// Load and pre-process gaussian data from .bin file
			await octreeGeometryLoader.load(node, octreeFileUrl);

			for (let cid = 0; cid < 8; cid++) {
				const child = node.children[cid];
				if (child) {
					// push child node into depth queue
					queue.push({ node: child, level: level + 1 });
				}
			}
		}
		
	}

	await new Promise((resolve) => setTimeout(resolve, 200));
	// init gaussian data and send to worker
	// 6*4 + 4 + 4 = 8*4
	// XYZ - Position (Float32)
	// XYZ - Scale (Float32)
	// RGBA - colors (uint8)
	// IJKL - quaternion/rot (uint8)
	gaussianSplats.rowLength = 3 * 4 + 3 * 4 + 4 + 4
	gaussianSplats.baseBuffer = await loadOctreeGeometry(octreeGeometry.root)
	gaussianSplats.baseVertexCount = gaussianSplats.baseBuffer.byteLength / gaussianSplats.rowLength
	// Send gaussian data to the worker
	worker.postMessage({
		buffer: gaussianSplats.baseBuffer,
		vertexCount: gaussianSplats.baseVertexCount,
	})

}
var progressTextDom = document.getElementById("progress-text");
var progressBarDom = document.getElementById("progressBar");
const chunkSize = 0.05 * 1024 * 1024 * 1024;//1G

const startLoadScene = (blob) => {
	octreeFileUrl = URL.createObjectURL(blob);
	progressBarDom.style.width = `100%`;
	loadScene()

	URL.revokeObjectURL(blob);
}
async function getFileData() {
	let data = await fetchSetFile(dataSource.octreeUrl)
	const mergedBlob = await mergeChunks(data);
	startLoadScene(mergedBlob)
	
}
async function fetchSetFile(url) {
	const totalSize = await getTotalSize(url);
	const chunkCount = Math.ceil(totalSize / chunkSize);
	let chunks = []
	let totalLoaded = 0; // Track the total bytes loaded

	for (let i = 0; i < chunkCount; i++) {
		const start = i * chunkSize;
		const end = Math.min(totalSize, (i + 1) * chunkSize) - 1;

		const xhr = new XMLHttpRequest();
		xhr.open('GET', url);
		xhr.setRequestHeader('Range', `bytes=${start}-${end}`);
		xhr.responseType = 'arraybuffer';

		const promise = new Promise((resolve, reject) => {
			xhr.onload = async () => {
				// 当从客户端发送Range范围标头以只请求资源的一部分时，将使用此响应代码。
				// 此处为读取一个chunk size的数据
				if (xhr.status === 206) {
					let blob = new Blob([xhr.response]);
					let url = URL.createObjectURL(blob);
					// await storeChunkInIndexedDB(blob, i);
					resolve(blob);
				} else {
					reject(new Error(`Failed to fetch chunk: ${xhr.statusText}`));
				}
			};
			xhr.onerror = () => {
				reject(new Error('Failed to fetch chunk: Network error'));
			};
		});
		// 
		xhr.addEventListener('progress', function (event) {
			if (event.lengthComputable) {
				// const progress = Math.round((totalLoaded / event.total) * 100);
				const progress = Math.round((i / chunkCount) * 100);
				totalLoaded = Math.max(totalLoaded, progress);
				progressBarDom.style.width = `${totalLoaded}%`;
				// progressTextDom.innerHTML = `Load Chunk ${i} Resources ${(dataSource.size * progress / 100).toFixed(2)}MB`;
				progressTextDom.innerHTML = `Load Resources ${(dataSource.size * totalLoaded / 100).toFixed(2)}MB`;
				if (progress === 100) {
					progressTextDom.innerText = `File Loaded Successfully, Parsing...`;
				}
			} else {
				console.log('Download progress: Not computable');
			}
		});

		xhr.send();
		chunks.push(promise);
	}

	return Promise.all(chunks);
}

async function getTotalSize(url) {
	const xhr = new XMLHttpRequest();
	xhr.open('HEAD', url);
	xhr.send();
	await new Promise(resolve => xhr.onload = resolve);
	return parseInt(xhr.getResponseHeader('Content-Length'));
}

async function storeChunkInIndexedDB(chunk, index) {
	const dbName = `${dataSource.name}_${index}`;
	await localforage.setItem(dbName, chunk);
}

// 获取分片数据并合并
async function mergeChunks(data) {
	const chunks = [];
	let index = 0;
	while (true) {
		let chunk
		if (data) {
			chunk = data[index]
		}
		else {
			chunk = await localforage.getItem(`${dataSource.name}_${index}`);
		}
		if (!chunk) break;
		chunks.push(chunk);
		index++;
	}
	const blob = new Blob(chunks);
	return blob;
}

getFileData()
// fetchSetFile(dataSource.octreeUrl)

async function updateGaussianDefault() {
	worker.postMessage({
		buffer: gaussianSplats.baseBuffer,
		vertexCount: gaussianSplats.baseVertexCount,
	})
}

function mergeBuffer(baseBuffer, extraBuffer) {
	// 1. 创建一个新的缓冲区，大小为 baseBuffer 和 extraBuffer 的总和
	let mergedBuffer = new ArrayBuffer(baseBuffer.byteLength + extraBuffer.byteLength);
	const mergedView = new Uint8Array(mergedBuffer);

    mergedView.set(new Uint8Array(baseBuffer), 0);
    mergedView.set(new Uint8Array(extraBuffer), baseBuffer.byteLength);

	return mergedBuffer;
}

async function updateGaussianByView(viewMatrix, projectionMatrix, maxLevel, maxCount) {
	// console.log(update_count, maxLevel);
	update_count += 1;
	if (update_count > 1) return;

	const start = performance.now();

	stopReading = false;
	gaussianSplats.extraVertexCount = 0;

	updateGaussianDefault();
	if (maxLevel === settings.baseLevel) { 
		update_count = 0;
		worker.postMessage({
			buffer: gaussianSplats.baseBuffer,
			vertexCount: gaussianSplats.baseVertexCount,
		})
		return; 
	}

	let ZDepthMax = settings.depthMax;
	let nodeList = [];
	let isOverMaxLimit = false;
	const maxConcurrentLoads = 16; // 最大并行加载数
	const taskQueue = []; // 存储待处理任务的队列
	
	for (let base_index = 0; base_index < baseLevelQueue.length; base_index++) {
		let queue = [];
	
		let { ZDepth, visibility } = markCubeVisibility(viewMatrix, projectionMatrix, baseLevelQueue[base_index].node);
		queue.push({ node: baseLevelQueue[base_index].node, level: baseLevelQueue[base_index].level, ZDepth: ZDepth, visibility: visibility });
		
		if (!visibility) {
			queue.pop();
		}

		if (queue.length > 0 && isOverMaxLimit) {
			const { node, level, ZDepth, visibility } = queue.shift();
			node.visibility = visibility;
			taskQueue.push(async () => {
				if (!node.geometry) {
					await octreeGeometryLoader.load(node, octreeFileUrl);
					// console.log(`load ${node.numPoints} completed!`) // 加载点云
				}
			});
			gaussianSplats.extraVertexCount += node.numPoints;
			nodeList.push(node);
		}
		
		while (queue.length > 0 && !isOverMaxLimit) {
			const { node, level, ZDepth, visibility } = queue.shift();
			node.reading = false;

			if (level <= maxLevel) {
				const beta = Math.log(maxLevel + 1);
				const actLevel = Math.floor(Math.min(Math.max((maxLevel + 1) * Math.exp(-1.0 * beta * Math.abs(ZDepth) / ZDepthMax), settings.baseLevel), maxLevel));
				node.visibility = level >= actLevel ? visibility : false;
	
				if (node.visibility && !node.reading) {
					// 将加载任务加入任务队列
					taskQueue.push(async () => {
						if (!node.geometry) {
							await octreeGeometryLoader.load(node, octreeFileUrl);
							// console.log(`load ${node.numPoints} completed!`) // 加载点云
						}
					});

					gaussianSplats.extraVertexCount += node.numPoints;
					nodeList.push(node);

					// 检查是否达到加载限制
					if (gaussianSplats.extraVertexCount > maxCount) {
						isOverMaxLimit = true;
						break;
					}
				} else if (level < maxLevel) {
				// 添加子节点到队列		
					for (let cid = 0; cid < 8; cid++) {
						const child = node.children[cid];
						if (child) {
							let { ZDepth, visibility } = markCubeVisibility(viewMatrix, projectionMatrix, child);
							if (visibility) {
								queue.push({ node: child, level: level + 1, ZDepth: ZDepth, visibility: visibility });
							}
						}
					}
				}
			}
		}
	}
	
	// 批量执行任务队列中的任务，控制并行加载的数量
	async function executeTaskQueue(taskQueue, maxConcurrentLoads) {
		const taskPromises = [];
		for (const task of taskQueue) {
			const taskPromise = task().catch(err => {
				console.error("Error loading node geometry:", err);
			});
			taskPromises.push(taskPromise);
	
			if (taskPromises.length >= maxConcurrentLoads) {
				// 等待当前并行任务完成后再启动新的任务
				await Promise.race(taskPromises);
				taskPromises.splice(taskPromises.findIndex(p => p.isResolved), 1); // 移除已完成的任务
			}
		}
		// 等待所有任务完成
		await Promise.all(taskPromises);
	}
	
	// 执行任务队列
	await executeTaskQueue(taskQueue, maxConcurrentLoads);

	const end_count = performance.now();
	const count_time = `${((end_count - start) / 1000).toFixed(3)}s`
	// Send gaussian data to the worker
	gaussianSplats.extraBuffer = new ArrayBuffer(gaussianSplats.rowLength * gaussianSplats.extraVertexCount);
	gaussianSplats.loadedCount = 0;
	gaussianSplats.lastloadedCount = 0;

	let campos = [viewMatrix[2], viewMatrix[6], viewMatrix[10]];
	// await readGaussianFromNode(octreeGeometry.root, gaussianSplats, campos, 0); // put something to extrabuffer
	await readGaussianFromNodeListMultiThread(nodeList, gaussianSplats, campos); // put something to extrabuffer
	const end_load = performance.now()
	const load_time = `${((end_load - end_count) / 1000).toFixed(3)}s`

	let merge_Buffer = mergeBuffer(gaussianSplats.baseBuffer, gaussianSplats.extraBuffer);
	
	// console.log(merge_Buffer); // true
	let merge_count = gaussianSplats.extraVertexCount + gaussianSplats.baseVertexCount;
	worker.postMessage({
		buffer: merge_Buffer,
		vertexCount: merge_count,
	})
	const total_time = `${((performance.now() - start) / 1000).toFixed(3)}s`
	progressTextDom.innerHTML = ``;
	reloadLod = false;
	console.log(`[Loader] count ${gaussianSplats.extraVertexCount} gaussians in ${count_time}.`)
	console.log(`[Loader] load ${gaussianSplats.extraVertexCount} gaussians in ${load_time}.`)
	console.log(`[Loader] total in ${total_time}.`)
	update_count = 0;
	gaussianSplats.extraVertexCount = 0;
	gaussianSplats.extraBuffer = null;

}

// async function updateGaussianByView(viewMatrix, projectionMatrix, maxLevel, maxCount) {
// 	// console.log(update_count, maxLevel);
// 	update_count += 1;
// 	if (update_count > 1) return;

// 	const start = performance.now();

// 	stopReading = false;
// 	gaussianSplats.extraVertexCount = 0;

// 	updateGaussianDefault();
// 	if (maxLevel === settings.baseLevel) { 
// 		update_count = 0;
// 		worker.postMessage({
// 			buffer: gaussianSplats.baseBuffer,
// 			vertexCount: gaussianSplats.baseVertexCount,
// 		})
// 		return; 
// 	}

// 	let ZDepthMax = settings.depthMax;
// 	let nodeList = [];
// 	let isOverMaxLimit = false;
	
// 	for (let base_index = 0; base_index < baseLevelQueue.length; base_index++) {
// 		let queue = [];
	
// 		let { ZDepth, visibility } = markCubeVisibility(viewMatrix, projectionMatrix, baseLevelQueue[base_index].node);
// 		queue.push({ node: baseLevelQueue[base_index].node, level: baseLevelQueue[base_index].level, ZDepth: ZDepth, visibility: visibility });
		
// 		if (!visibility) {
// 			queue.pop();
// 		}

// 		if (queue.length > 0 && isOverMaxLimit) {
// 			const { node, level, ZDepth, visibility } = queue.shift();
// 			node.visibility = visibility;
// 			if (!node.geometry) {
// 				await octreeGeometryLoader.load(node, octreeFileUrl);
// 				// console.log(`load ${node.numPoints} completed!`) // 加载点云
// 			}

// 			gaussianSplats.extraVertexCount += node.numPoints;
// 			nodeList.push(node);
// 		}
		
// 		while (queue.length > 0 && !isOverMaxLimit) {
// 			const { node, level, ZDepth, visibility } = queue.shift();
// 			node.reading = false;

// 			if (level <= maxLevel) {
// 				const beta = Math.log(maxLevel + 1);
// 				const actLevel = Math.floor(Math.min(Math.max((maxLevel + 1) * Math.exp(-1.0 * beta * Math.abs(ZDepth) / ZDepthMax), settings.baseLevel), maxLevel));
// 				node.visibility = level >= actLevel ? visibility : false;
	
// 				if (node.visibility && !node.reading) {
// 					// 将加载任务加入任务队列
// 					if (!node.geometry) {
// 						await octreeGeometryLoader.load(node, octreeFileUrl);
// 						// console.log(`load ${node.numPoints} completed!`) // 加载点云
// 					}

// 					gaussianSplats.extraVertexCount += node.numPoints;
// 					nodeList.push(node);

// 					// 检查是否达到加载限制
// 					if (gaussianSplats.extraVertexCount > maxCount) {
// 						isOverMaxLimit = true;
// 						break;
// 					}
// 				} else if (level < maxLevel) {
// 				// 添加子节点到队列		
// 					for (let cid = 0; cid < 8; cid++) {
// 						const child = node.children[cid];
// 						if (child) {
// 							let { ZDepth, visibility } = markCubeVisibility(viewMatrix, projectionMatrix, child);
// 							if (visibility) {
// 								queue.push({ node: child, level: level + 1, ZDepth: ZDepth, visibility: visibility });
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
	
// 	const end_count = performance.now();
// 	const count_time = `${((end_count - start) / 1000).toFixed(3)}s`
// 	// Send gaussian data to the worker
// 	gaussianSplats.extraBuffer = new ArrayBuffer(gaussianSplats.rowLength * gaussianSplats.extraVertexCount);
// 	gaussianSplats.loadedCount = 0;
// 	gaussianSplats.lastloadedCount = 0;

// 	let campos = [viewMatrix[2], viewMatrix[6], viewMatrix[10]];
// 	// await readGaussianFromNode(octreeGeometry.root, gaussianSplats, campos, 0); // put something to extrabuffer
// 	await readGaussianFromNodeList(nodeList, gaussianSplats, campos); // put something to extrabuffer
// 	const end_load = performance.now()
// 	const load_time = `${((end_load - end_count) / 1000).toFixed(3)}s`

// 	let merge_Buffer = mergeBuffer(gaussianSplats.baseBuffer, gaussianSplats.extraBuffer);
	
// 	// console.log(merge_Buffer); // true
// 	let merge_count = gaussianSplats.extraVertexCount + gaussianSplats.baseVertexCount;
// 	worker.postMessage({
// 		buffer: merge_Buffer,
// 		vertexCount: merge_count,
// 	})
// 	const total_time = `${((performance.now() - start) / 1000).toFixed(3)}s`
// 	progressTextDom.innerHTML = ``;
// 	reloadLod = false;
// 	console.log(`[Loader] count ${gaussianSplats.extraVertexCount} gaussians in ${count_time}.`)
// 	console.log(`[Loader] load ${gaussianSplats.extraVertexCount} gaussians in ${load_time}.`)
// 	console.log(`[Loader] total in ${total_time}.`)
// 	update_count = 0;
// 	gaussianSplats.extraVertexCount = 0;
// 	gaussianSplats.extraBuffer = null;

// }

async function readGaussianFromNode(node, gaussianSplats, campos, level) {
	if (node.visibility && node.loaded && !node.reading) {
		//set reading flag
		node.reading = true;
		const currentCount = node.numPoints
		for (let i = 0; i < currentCount; i++) {
			// check buffer size
			if (!isBufferLargeEnough(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3, 4)) return;
			// read into buffer
			const positions = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3);
			const scales = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3, 3);
			const rgbas = new Uint8ClampedArray(
				gaussianSplats.extraBuffer,
				gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3,
				4,
			);
			const rots = new Uint8ClampedArray(
				gaussianSplats.extraBuffer,
				gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3 + 4,
				4,
			);

			let { position, harmonic, opacity, scale, rotation } = nodeView(i, node)
			// Normalize quaternion
			let length2 = 0

			for (let j = 0; j < 4; j++)
				length2 += rotation[j] * rotation[j]

			const length = Math.sqrt(length2)

			rotation = rotation.map(v => (v / length) * 128 + 128)

			rots[0] = rotation[0];
			rots[1] = rotation[1];
			rots[2] = rotation[2];
			rots[3] = rotation[3];

			// Exponentiate scale
			scale = scale.map(v => Math.exp(v))

			scales[0] = scale[0];
			scales[1] = scale[1];
			scales[2] = scale[2];

			// const SH_C0 = 0.28209479177387814
			// rgbas[0] = (0.5 + SH_C0 * harmonic[0]) * 255;
			// rgbas[1] = (0.5 + SH_C0 * harmonic[1]) * 255;
			// rgbas[2] = (0.5 + SH_C0 * harmonic[2]) * 255;

			let color = computeColorFromSH(settings.shDegree, position, campos, harmonic)
			rgbas[0] = color.x * 255;
			rgbas[1] = color.y * 255;
			rgbas[2] = color.z * 255;

			// Activate alpha
			const sigmoid = (m1) => 1 / (1 + Math.exp(-m1))
			rgbas[3] = sigmoid(opacity) * 255;

			positions[0] = position[0];
			positions[1] = position[1];
			positions[2] = position[2];

			gaussianSplats.loadedCount++;
		}
	}
	// set the children of the current node
	for (let cid = 0; cid < 8; cid++) {
		const child = node.children[cid];
		if (child) {
			await readGaussianFromNode(child, gaussianSplats, campos, level + 1)
		}
	}
}

// read gaussian data from octree node List saved
async function readGaussianFromNodeList(nodeList, gaussianSplats, campos) {
	gaussianSplats.loadedCount = 0;
	for (let idx = 0; idx < nodeList.length; idx++) {
		if (nodeList[idx].visibility && nodeList[idx].loaded && !nodeList[idx].reading) {
			//set reading flag
			nodeList[idx].reading = true;
			const currentCount = nodeList[idx].numPoints
			for (let i = 0; i < currentCount; i++) {
				// check buffer size
				if (!isBufferLargeEnough(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3, 4)) return;
				// read into buffer
				const positions = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3);
				const scales = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3, 3);
				const rgbas = new Uint8ClampedArray(
					gaussianSplats.extraBuffer,
					gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3,
					4,
				);
				const rots = new Uint8ClampedArray(
					gaussianSplats.extraBuffer,
					gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3 + 4,
					4,
				);

				let { position, harmonic, opacity, scale, rotation } = nodeView(i, nodeList[idx])
				// Normalize quaternion
				let length2 = 0

				for (let j = 0; j < 4; j++)
					length2 += rotation[j] * rotation[j]

				const length = Math.sqrt(length2)

				rotation = rotation.map(v => (v / length) * 128 + 128)

				rots[0] = rotation[0];
				rots[1] = rotation[1];
				rots[2] = rotation[2];
				rots[3] = rotation[3];

				// Exponentiate scale
				scale = scale.map(v => Math.exp(v))

				scales[0] = scale[0];
				scales[1] = scale[1];
				scales[2] = scale[2];

				// const SH_C0 = 0.28209479177387814
				// rgbas[0] = (0.5 + SH_C0 * harmonic[0]) * 255;
				// rgbas[1] = (0.5 + SH_C0 * harmonic[1]) * 255;
				// rgbas[2] = (0.5 + SH_C0 * harmonic[2]) * 255;

				let color = computeColorFromSH(settings.shDegree, position, campos, harmonic)
				rgbas[0] = color.x * 255;
				rgbas[1] = color.y * 255;
				rgbas[2] = color.z * 255;

				// Activate alpha
				const sigmoid = (m1) => 1 / (1 + Math.exp(-m1))
				rgbas[3] = sigmoid(opacity) * 255;

				positions[0] = position[0];
				positions[1] = position[1];
				positions[2] = position[2];

				gaussianSplats.loadedCount++;
			}
		}
	}
}
async function readGaussianFromNodeListMultiThread(nodeList, gaussianSplats, campos) {
    gaussianSplats.loadedCount = 0;

    // 处理每个节点的异步任务
    const processNode = async (node) => {
        if (node.visibility && node.loaded && !node.reading) {
            node.reading = true;
            const currentCount = node.numPoints;

            for (let i = 0; i < currentCount; i++) {
                // 检查缓冲区大小
                if (!isBufferLargeEnough(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3, 4)) return;

                // 获取当前点的数据
                const positions = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength, 3);
                const scales = new Float32Array(gaussianSplats.extraBuffer, gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3, 3);
                const rgbas = new Uint8ClampedArray(
                    gaussianSplats.extraBuffer,
                    gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3,
                    4,
                );
                const rots = new Uint8ClampedArray(
                    gaussianSplats.extraBuffer,
                    gaussianSplats.loadedCount * gaussianSplats.rowLength + 4 * 3 + 4 * 3 + 4,
                    4,
                );

                // 获取点的数据
                let { position, harmonic, opacity, scale, rotation } = nodeView(i, node);

                // 归一化四元数
                let length2 = 0;
                for (let j = 0; j < 4; j++) length2 += rotation[j] * rotation[j];
                const length = Math.sqrt(length2);
                rotation = rotation.map(v => (v / length) * 128 + 128);

                // 将归一化后的四元数写入 rots
                rots[0] = rotation[0];
                rots[1] = rotation[1];
                rots[2] = rotation[2];
                rots[3] = rotation[3];

                // 指数化缩放因子
                scale = scale.map(v => Math.exp(v));

                // 将缩放因子写入 scales
                scales[0] = scale[0];
                scales[1] = scale[1];
                scales[2] = scale[2];

                // 计算颜色
                let color = computeColorFromSH(settings.shDegree, position, campos, harmonic);
                rgbas[0] = color.x * 255;
                rgbas[1] = color.y * 255;
                rgbas[2] = color.z * 255;

                // 使用 Sigmoid 函数计算透明度
                const sigmoid = (m1) => 1 / (1 + Math.exp(-m1));
                rgbas[3] = sigmoid(opacity) * 255;

                // 写入位置
                positions[0] = position[0];
                positions[1] = position[1];
                positions[2] = position[2];

                // 更新计数器
                gaussianSplats.loadedCount++;
            }
        }
    };

    // 创建一个数组来保存所有的异步任务
    const processTasks = nodeList.map(processNode);

    // 等待所有的任务完成
    await Promise.all(processTasks);
}


function markCubeVisibility(viewMatrix, projectionMatrix, node) {
	// get 8 boundary points of the cube
	const pmin = node.boundingBox.min;
	const pmax = node.boundingBox.max;
	const p000 = [pmin.x, pmin.y, pmin.z];
	const p001 = [pmin.x, pmin.y, pmax.z];
	const p010 = [pmin.x, pmax.y, pmin.z];
	const p011 = [pmin.x, pmax.y, pmax.z];
	const p100 = [pmax.x, pmin.y, pmin.z];
	const p101 = [pmax.x, pmin.y, pmax.z];
	const p110 = [pmax.x, pmax.y, pmin.z];
	const p111 = [pmax.x, pmax.y, pmax.z];
	const points = [p000, p001, p010, p011, p100, p101, p110, p111];


	const center = [node.boundingSphere.center.x, node.boundingSphere.center.y, node.boundingSphere.center.z];
	const { ZDepth, visibility } = markPointVisibility(viewMatrix, projectionMatrix, center);

	if  (visibility) return { ZDepth, visibility };

	// for (let i = 0; i < 8; i++) {
	// 	const { ZDepth, visibility } = markPointVisibility(viewMatrix, projectionMatrix, points[i]);
	// 	if (visibility) return { ZDepth, visibility };
	// }

	return { ZDepth, visibility };
}

function markPointVisibility(viewMatrix, projectionMatrix, centerPos) {
	const x = centerPos[0];
	const y = centerPos[1];
	const z = centerPos[2];

	// Convert to camera space
	const ZDepth = x * viewMatrix[2] + y * viewMatrix[6] + z * viewMatrix[10] + viewMatrix[14];

	// Get View Projection Matrix
	const viewProjMatrix = multiply4(projectionMatrix, viewMatrix);

	// Convert to clip space
	let clipX = x * viewProjMatrix[0] + y * viewProjMatrix[4] + z * viewProjMatrix[8] + viewProjMatrix[12];
	let clipY = x * viewProjMatrix[1] + y * viewProjMatrix[5] + z * viewProjMatrix[9] + viewProjMatrix[13];
	let clipZ = x * viewProjMatrix[2] + y * viewProjMatrix[6] + z * viewProjMatrix[10] + viewProjMatrix[14];
	let clipW = x * viewProjMatrix[3] + y * viewProjMatrix[7] + z * viewProjMatrix[11] + viewProjMatrix[15];

	clipX = clipX / (clipW + 0.000001);
	clipY = clipY / (clipW + 0.000001);
	clipZ = clipZ / (clipW + 0.000001);

	let visibility = false;
	// Check if the point is within the camera's field of view
	// modified the visible depth to 3
	if (ZDepth > -3 && clipX > -1.5 && clipX < 1.5 && clipY > -1.5 && clipY < 1.5) {
		visibility = true;
	}

	return { ZDepth, visibility };
}

function isBufferLargeEnough(buffer, offset, length, byteSize) {
	let requiredBytes = length * byteSize;
	return buffer.byteLength >= offset + requiredBytes;
}

function computeColorFromSH(deg, position, campos, harmonic) {
	const SH_C0 = 0.28209479177387814;
	const SH_C1 = 0.4886025119029199;
	const SH_C2 = [1.0925484305920792,
		-1.0925484305920792,
		0.31539156525252005,
		-1.0925484305920792,
		0.5462742152960396];

	const SH_C3 = [
		-0.5900435899266435,
		2.890611442640554,
		-0.4570457994644658,
		0.3731763325901154,
		-0.4570457994644658,
		1.445305721320277,
		-0.5900435899266435
	];

	// Helper function to normalize a vector
	function normalize(vec) {
		let len = Math.sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
		return { x: vec.x / len, y: vec.y / len, z: vec.z / len };
	}

	// Helper function to compute the maximum of a vector and 0
	function maxVec(vec) {
		return { x: Math.max(vec.x, 0), y: Math.max(vec.y, 0), z: Math.max(vec.z, 0) };
	}

	// let pos = { x: position[0], y: position[1], z: position[2] };
	let dir = { x: position[0] - campos[0], y: position[1] - campos[1], z: position[2] - campos[2] };
	dir = normalize(dir);
	// dir = { x: -dir.x, y: -dir.y, z: -dir.z };

	let sh = harmonic;

	let result = { x: SH_C0 * sh[0], y: SH_C0 * sh[1], z: SH_C0 * sh[2] };

	if (deg > 0) {
		let x = dir.x;
		let y = dir.y;
		let z = dir.z;
		result = {
			x: result.x - SH_C1 * y * sh[3] + SH_C1 * z * sh[6] - SH_C1 * x * sh[9],
			y: result.y - SH_C1 * y * sh[4] + SH_C1 * z * sh[7] - SH_C1 * x * sh[10],
			z: result.z - SH_C1 * y * sh[5] + SH_C1 * z * sh[8] - SH_C1 * x * sh[11]
		};
		// console.log(dir, sh)

		if (deg > 1) {
			let xx = x * x, yy = y * y, zz = z * z;
			let xy = x * y, yz = y * z, xz = x * z;
			result = {
				x: result.x + SH_C2[0] * xy * sh[12] + SH_C2[1] * yz * sh[15] + SH_C2[2] * (2.0 * zz - xx - yy) * sh[18] + SH_C2[3] * xz * sh[21] + SH_C2[4] * (xx - yy) * sh[24],
				y: result.y + SH_C2[0] * xy * sh[13] + SH_C2[1] * yz * sh[16] + SH_C2[2] * (2.0 * zz - xx - yy) * sh[19] + SH_C2[3] * xz * sh[22] + SH_C2[4] * (xx - yy) * sh[25],
				z: result.z + SH_C2[0] * xy * sh[14] + SH_C2[1] * yz * sh[17] + SH_C2[2] * (2.0 * zz - xx - yy) * sh[20] + SH_C2[3] * xz * sh[23] + SH_C2[4] * (xx - yy) * sh[26]
			};

			if (deg > 2) {
				// Add the third-degree SH terms
				// Similar to the second-degree terms, add the necessary calculations
			}
		}
	}
	result.x += 0.5;
	result.y += 0.5;
	result.z += 0.5;

	return maxVec(result);
	// return result;
}
