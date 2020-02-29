var size;
var canvas;
var ctx;
var shrink = .75;
var padding = .025;
var paint_interval = 30; // ms
var paint_timer;
var rad = .02;
var nmols = 20;
var kT = 10;
var vrms = Math.sqrt(kT*1E-6); // pixels per second
var mols = [];
var energy = 0.;
var efield = .001;
var drag = 0.01;
var runningaverage = 200;
var avgmoment = null; 
for (var i = 0; i < nmols; i++) {
    x = 2.*rad + (1-2.*rad) * Math.random();
    y = 2.*rad + (1-2.*rad) * Math.random();
    theta = 2.*Math.PI * (Math.random() - .5);
    vx = Math.random() - .5;
    vy = Math.random() - .5;
    omega = 1./Math.sqrt(2)/rad*(Math.random() - .5);
    energy += vx**2 + vy**2 + (rad*omega)**2;
    mols.push(
	[x,y,theta,vx,vy,omega]
    );
}
energy_scale = nmols * vrms**2 / energy;
velocity_scale = Math.sqrt(energy_scale)
for (var mol of mols) {
    for (var j = 0; j < 3; j++) {
	mol[3+j] *= velocity_scale;
    }
}

function draw_circle(x,y,rad,strokecolor,fillcolor) {
    ctx.beginPath();
    ctx.arc(scalex(x),scaley(y),scale(rad),0.,2.*Math.PI);
    ctx.fillStyle = fillcolor;
    ctx.fill();
    ctx.strokeStyle = strokecolor;
    ctx.lineWidth = scale(rad)/5.;
    ctx.stroke();
}

function draw_plus(x,y,width,strokecolor) {
    ctx.beginPath();
    px = scalex(x);
    py = scaley(y);
    pw = scale(width);
    ctx.moveTo(px-pw/2,py);
    ctx.lineTo(px+pw/2,py);
    ctx.moveTo(px,py-pw/2);
    ctx.lineTo(px,py+pw/2);
    ctx.strokeStyle = strokecolor;
    ctx.lineWidth = pw / 4.;
    ctx.stroke();
}

function draw_minus(x,y,width,strokecolor) {
    ctx.beginPath();
    px = scalex(x);
    py = scaley(y);
    pw = scale(width);
    ctx.moveTo(px-pw/2,py);
    ctx.lineTo(px+pw/2,py);
    ctx.strokeStyle = strokecolor;
    ctx.lineWidth = pw / 4.;
    ctx.stroke();
}

function draw_mol(x,y,theta) {
    var rct = rad*Math.cos(theta);
    var rst = rad*Math.sin(theta);
    draw_circle(
	x + rct,
	y + rst,
	rad,
	'#444444',
	'#CC8888',
    );
    draw_plus(x + rct,y + rst,.8*rad,'#444444');
    draw_circle(
	x - rct,
	y - rst,
	rad,
	'#444444',
	'#8888CC',
    );
    draw_minus(x - rct,y - rst,.8*rad,'#444444');
}

function draw_box(xll,yll,xur,yur,strokecolor,fillcolor) {
    ctx.beginPath();
    ctx.moveTo(scalex(xll),scaley(yll));
    ctx.lineTo(scalex(xur),scaley(yll));
    ctx.lineTo(scalex(xur),scaley(yur));
    ctx.lineTo(scalex(xll),scaley(yur));
    ctx.closePath();
    ctx.fillStyle = fillcolor;
    ctx.fill();
    ctx.strokeStyle = strokecolor;    
    ctx.stroke();
}

function scale(x) {
    return size*x;
}

function scalex(x) {
    return size*(x+padding);
}

function scaley(y) {
    return size*(1+padding-y);
}

function left_wall_plus(x,y,theta) {
    return x - rad * (1 - Math.cos(theta)) < 0
}

function left_wall_minus(x,y,theta) {
    return x - rad * (1 + Math.cos(theta)) < 0
}

function right_wall_plus(x,y,theta) {
    return x + rad * (1 + Math.cos(theta)) > 1
}

function right_wall_minus(x,y,theta) {
    return x + rad * (1 - Math.cos(theta)) > 1
}

function bottom_wall_plus(x,y,theta) {
    return y - rad * (1 - Math.sin(theta)) < 0
}

function bottom_wall_minus(x,y,theta) {
    return y - rad * (1 + Math.sin(theta)) < 0
}

function top_wall_plus(x,y,theta) {
    return y + rad * (1 + Math.sin(theta)) > 1
}

function top_wall_minus(x,y,theta) {
    return y + rad * (1 - Math.sin(theta)) > 1
}

function check_collision(x,y,theta) {
    for (
	var cond of [
	    left_wall_plus,
	    left_wall_minus,
	    right_wall_plus,
	    right_wall_minus,
	    bottom_wall_plus,
	    bottom_wall_minus,
	    top_wall_plus,
	    top_wall_minus	   
	]
    ) {
	if (cond(x,y,theta)) {
	    return true;
	}
    }
    return false;
}

function vrx(mol,dx,dy) {
    var omega = mol[5];
    // omega zhat X ( dx xhat + dy yhat ) dot xhat
    return -omega * dy
}

function vry(omega,dx,dy) {
    // omega zhat X ( dx xhat + dy yhat ) dot yhat
    return omega * dx
}


function draw() {
    canvas.height = size*(1.+2.*padding);
    canvas.width = size*(1.+2.*padding);    
    draw_box(0.,0.,1.,1.,'black','white');
    draw_box(0.,-padding,1.,0.,'black',	'#CC8888');
    draw_box(0.,1.,1.,1.+padding,'black','#8888CC');
    rot_energy = 0;
    trans_energy = 0;
    moment = 0;
    for (var i = 0; i < nmols; i++) {
	var mol = mols[i];
	var x = mol[0];
	var y = mol[1];
	var theta = mol[2];
	var sintheta = Math.sin(theta);
	var costheta = Math.cos(theta);
	var vx = mol[3];
	var vy = mol[4];
	var omega = mol[5];
	draw_mol(x,y,theta);
	// left side
	if (x-rad*(1+Math.abs(Math.cos(theta))) < 0) {	    
	    var st = sintheta;
	    st *= Math.abs(theta) > Math.PI / 2 ? +1 : -1;
	    var vrelminus = vx - omega*rad*st;
	    if (vrelminus < 0) {
		var j = -2. * vrelminus / (1. + st**2);
		mol[3] += j
		mol[5] += -st*j/rad;
	    }
	}
	// right side
	if (x+rad*(1+Math.abs(Math.cos(theta))) > 1.) {
	    var st = sintheta;
	    st *= Math.abs(theta) < Math.PI / 2 ? +1 : -1;
	    var vrelminus = vx - omega*rad*st;
	    if (vrelminus > 0) {
		var j = -2. * vrelminus / (1. + st**2);
		mol[3] += j
		mol[5] += -st*j/rad;
	    }
	}
	// bottom side
	if (y-rad*(1+Math.abs(Math.sin(theta))) < 0.) {
	    var t = theta < 0 ? theta : theta + Math.PI;
	    var ct = Math.cos(t);
	    var vrelminus = vy + omega*rad*ct;
	    if (vrelminus < 0) {
		var j = -2. * vrelminus / (1. + ct**2);
		mol[4] += j
		mol[5] += ct*j/rad;
	    }
	}
	// top side
	if (y+rad*(1+Math.abs(Math.sin(theta))) > 1.) {	    
	    var t = theta > 0 ? theta : theta + Math.PI;
	    var ct = Math.cos(t);
	    var vrelminus = -vy - omega*rad*ct;
	    if (vrelminus < 0) {
		var j = -2. * vrelminus / (1. + ct**2);
		mol[4] += -j
		mol[5] += -ct*j/rad;
	    }
	}
	// collisions
	var xp = x + rad*costheta;
	var yp = y + rad*sintheta;
	var xm = x - rad*costheta;
	var ym = y - rad*sintheta;
	for (var jjj = 0; jjj < i; jjj++) {
	    var Mol = mols[jjj];
	    var X = Mol[0];
	    var Y = Mol[1];
	    var Theta = Mol[2];
	    var Sintheta = Math.sin(Theta);
	    var Costheta = Math.cos(Theta);
	    var Vx = Mol[3];
	    var Vy = Mol[4];
	    var Omega = Mol[5];
	    var Xp = X + rad*Costheta;
	    var Yp = Y + rad*Sintheta;
	    var Xm = X - rad*Costheta;
	    var Ym = Y - rad*Sintheta;
	    if (
		(xp-Xp)**2 + (yp-Yp)**2 < (2*rad)**2
		    ||
		    (xm-Xp)**2 + (ym-Yp)**2 < (2*rad)**2
		    ||
		    (xp-Xm)**2 + (yp-Ym)**2 < (2*rad)**2
		    ||
		    (xm-Xm)**2 + (ym-Ym)**2 < (2*rad)**2
	    ) {
		if (
		    (xp-Xp)**2 + (yp-Yp)**2 < (2*rad)**2
		) {
		    costheta *= +1;
		    sintheta *= +1;
		    Costheta *= +1;
		    Sintheta *= +1;
		}
		else if (
		    (xm-Xp)**2 + (ym-Yp)**2 < (2*rad)**2
		) {
		    costheta *= -1;
		    sintheta *= -1;
		    Costheta *= +1;
		    Sintheta *= +1;
		}
		else if (
		    (xp-Xm)**2 + (yp-Ym)**2 < (2*rad)**2
		) {
		    costheta *= +1;
		    sintheta *= +1;
		    Costheta *= -1;
		    Sintheta *= -1;
		}
		else if (
		    (xm-Xm)**2 + (ym-Ym)**2 < (2*rad)**2
		) {
		    costheta *= -1;
		    sintheta *= -1;
		    Costheta *= -1;
		    Sintheta *= -1;
		}
		
		// relative velocity of head centers
		var vhx = vx - Vx - rad * ( omega * sintheta - Omega * Sintheta);
		var vhy = vy - Vy + rad * ( omega * costheta - Omega * Costheta);
		var hx = xp - Xp;
		var hy = yp - Yp;
		if (hx * vhx + hy * vhy < 0) {
		    // separation decreasing, handle collision
		    var h = Math.sqrt(hx**2 + hy**2);
		    var nx = hx / h;
		    var ny = hy / h;
		    var vrelminus = vhx*nx + vhy*ny;
		    var rax = rad * costheta - hx / 2.;
		    var ray = rad * sintheta - hy / 2.;
		    var Rax = rad * Costheta + hx / 2.;
		    var Ray = rad * Sintheta + hy / 2.;
		    var j = -2.*vrelminus/(
			2. + 1./rad**2 * (
			    (rax**2 + ray**2) - (rax*nx + ray*ny)**2
			    + (Rax**2 + Ray**2) - (Rax*nx + Ray*ny)**2
			)
		    );
		    var jx = j * nx;
		    var jy = j * ny;
		    mol[3] += jx;
		    mol[4] += jy;
		    Mol[3] -= jx;
		    Mol[4] -= jy;
		    mol[5] += 1./rad**2*(rax*jy - ray*jx);
		    Mol[5] -= 1./rad**2*(Rax*jy - Ray*jx);
		}
	    }
	}
	mol[0] += mol[3];
	mol[1] += mol[4];
	mol[2] += mol[5];
	if (mol[2] > Math.PI) {
	    mol[2] -= 2* Math.PI;
	}
	else if (mol[2] < -Math.PI) {
	    mol[2] += 2* Math.PI;
	}
	mol[5] += efield*Math.cos(mol[2]) - drag*mol[5];
	trans_energy += mol[3]**2 + mol[4]**2;
	rot_energy += (rad*mol[5])**2;
	moment += Math.sin(mol[2]);
    }
    energy_scale = Math.max(
	(nmols * vrms**2 - rot_energy)/trans_energy,
	.1
    );    
    velocity_scale = Math.sqrt(energy_scale);
    var moment = moment / nmols;
    avgmoment = (
	avgmoment == null ? moment : (runningaverage-1) * avgmoment / runningaverage
    ) + moment / runningaverage;
    document.getElementById('alignment').innerText = Math.round(100 * avgmoment) + '%';
    for (var mol of mols) {
    	mol[3] *= velocity_scale;
    	mol[4] *= velocity_scale;
    }
}

function tick() {
    draw();
    setTimeout(tick,0);
}

function ontemperature(temperature) {
    var temperature = document.getElementById('temperature').value;
    kT = temperature;
    vrms = Math.sqrt(1e-6*kT);
    var energy = 0;
    for (mol of mols) {
	energy += mol[3]**2 + mol[4]**2 + (rad*mol[5])**2;
    }
    var energy_scale = nmols*vrms**2/energy;
    var velocity_scale = Math.sqrt(energy_scale);
    for (mol of mols) {
	for (var i = 0; i < 3; i++) {
	    mol[3+i] *= velocity_scale;
	}
    }
}

function onelectricfield(electricfield) {
    var electricfield = document.getElementById('electricfield').value;
    efield = electricfield*1e-4;
}

function on_load() {
    set_size();
    canvas = document.getElementById('canvas');
    ctx = canvas.getContext('2d');
    document.getElementById('temperature').onchange = ontemperature;
    document.getElementById('electricfield').onchange = onelectricfield;
    ontemperature();
    onelectricfield();
    tick();
};

function set_size() {
    size = shrink*Math.min(
	window.innerHeight,window.innerWidth
    );
}

window.onload = on_load;
window.onresize = set_size;
window.onkeypress = function (event) {on_key_press(event.key);};
