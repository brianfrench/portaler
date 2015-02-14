





frameRate(0);
mouseX=width>>1|7;
mouseY=height>>1|7;


var GF=function(){
    this.mul=function(a, x){
        for(var i=0, b=[]; i<a.length; i++){
            b[i]=a[i]*x;
        }
        return b;
    };
    this.add=function(a, b){
        for(var i=0, c=[]; i<a.length; i++){
            c[i]=a[i] + b[i];
        }
        return c;
    };
    this.sub=function(a, b){
        // return this.add(a, this.mul(b, -1));
        for(var i=0, c=[]; i<a.length; i++){
            c[i]=a[i] - b[i];
        }
        return c;
    };
    this.dot=function(a, b){
        for(var i=0, x=0; i<a.length; i++){
            x+=a[i]*b[i];
        }
        return x;
    };
    this.cross=function(a, b){
        for(var i=0, c=[]; i<3; i++){
            c[i]=a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3];
        }
        return c;
    };
    this.mag2=function(a){
        return this.dot(a, a);
    };
    this.norm=function(a){
        return this.mul(a, pow(this.mag2(a), -0.5));
    };
    this.lerp=function(a, b, t){
        return this.add(this.mul(a, 1 - t), this.mul(b, t));
    };
    
    this.pclip=function(poly, p, n){
        // if(n===[[0,0,0]]){return [];}
        if(this.mag2(n)<1e-9){return [];}
        var op=[];
        var m0, m1, v0=poly[poly.length - 1], v1;
        m0=this.dot(this.sub(v0, p), n);
        for(var i=0, j=0; i<poly.length; i++){
            var v1=poly[i].slice();
            m1=this.dot(this.sub(v1, p), n);
            if(m0*m1<0){
                op[j++]=this.lerp(v0, v1, m0/(m0 - m1));
            }
            if(m1>=0){
                op[j++]=v1.slice();
            }
            m0=m1;
            v0=v1;
        }
        return op;
    };
    this.cull=function(poly, o, beam){
        var op=poly.slice();
        for(var i=0; i<beam.length && op.length; i++){
            op=this.pclip(op, o, beam[i]);
        }
        return op;
    };
    this.beamify=function(poly, o){
        var beam=[], n=poly.length;
        if(n<3){return [[0,0,0]];}//[[-1,0,0],[1,0,0],];}
        var sgn=this.dot(
                    this.cross(this.sub(poly[1], poly[0]),
                               this.sub(poly[2], poly[1])),
                    this.sub(poly[0], o));
        sgn=sgn<0?1:-1;
        for(var i=0; i<n; i++){
            beam[i]=this.mul(this.cross(
                this.sub(poly[(i+1)%n], poly[i]),
                this.sub(poly[i], o)), sgn);
        }
        return beam;
    };
    
    // this.pih=function(p, h){
        
    // };
    
    this.qm=function(a, b){
        return [a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
                a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
                a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
                a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]];
    };
    this.qp=function(a){
        return [a[0], -a[1], -a[2], -a[3]];
    };
    this.qrot=function(q, p){
        var r=this.qm(this.qm(q, [0, p[0], p[1], p[2]]), this.qp(q));
        return [r[1], r[2], r[3]];
    };
};

var rq=[1, 0, 0, 0], vo=[1.1, 1.1, 1.1], vdo=[0.1, 0, 0], vhi=0;
rq=[2,1,1,1];
var mapon=1;



var kp=[], keyPressed=function(){
    kp[keyCode]=1;
    // if(keyCode===SHIFT){gtex.length=0;ntex=[0,6,0,0];}
    if(keyCode==='M'.charCodeAt(0)){mapon^=1;}
    
},keyReleased=function(){kp[keyCode]=0;};
for(var i=0; i<256; i++){kp[i]=0;}


var lms=0;
var phys=function(){
    var g=new GF();
    
    var ms=millis(), dt=lms ? (ms - lms)*1e-3 : 0;
    if(dt>0.3){dt=0.3;}
    lms=ms;
    
    // rq=g.norm(g.qm(rq,
    //     [1, -(mouseY - height/2)*1e-4, (mouseX - width/2)*1e-4, 0]));
    // vdo=g.mul(g.add(vdo, g.mul(g.qrot(rq, [0,0,1]), 
    //     (mouseIsPressed?0.01:0)*(mouseButton===RIGHT?-1:1))), 0.9);
    rq=g.norm(g.qm(rq,
        [1, (kp[UP] - kp[DOWN])*2e-2, (kp[RIGHT] - kp[LEFT])*2e-2, 0]));
    
    vdo=g.mul(vdo, pow(0.1, dt));
    vdo=g.add(vdo, g.mul(g.qrot(rq, [0,0,1]), 
        (kp['W'.charCodeAt(0)] - kp['S'.charCodeAt(0)])*10*dt));
    vdo=g.add(vdo, g.mul(g.qrot(rq, [1,0,0]), 
        (kp['D'.charCodeAt(0)] - kp['A'.charCodeAt(0)])*10*dt));
    vdo=g.add(vdo, g.mul(g.qrot(rq, [0,1,0]), 
        (kp['Q'.charCodeAt(0)] - kp['E'.charCodeAt(0)])*10*dt));
    
    vdo=g.add(vdo, [0,sin(ms*2e-1)*4e-3,0]);
    vo=g.add(vo, g.mul(vdo, dt));
};

var rscr=function(p,wid,hei, pt, po, h){
    var g=new this.GF();
    var vo=po[0], vf=po[1], vr=po[2], vd=po[3];
    var vu=pt[2], vv=pt[3], vn=pt[1], vt=pt[0];
    var ve=g.sub(po[0], vt);
    var siz=198, isiz=1/siz;
    // var h=~~0x495f3210;
    for(var y=-siz; y<siz; y++){
        var vb=g.add(vf, g.mul(vd, y*isiz));
        var b=g.dot(vn, vr)*isiz, a=g.dot(vn, vb) - siz*b;
        var d=(g.dot(ve, vu)*g.dot(vn, vr) + g.dot(vr, vu)*g.dot(vn, ve))*isiz, c=g.dot(ve, vu)*g.dot(vn, vb) + g.dot(vb, vu)*g.dot(vn, ve) - siz*d;
        var f=(g.dot(ve, vv)*g.dot(vn, vr) + g.dot(vr, vv)*g.dot(vn, ve))*isiz, e=g.dot(ve, vv)*g.dot(vn, vb) + g.dot(vb, vv)*g.dot(vn, ve) - siz*f;
        var l=80200 + y*wid - siz << 2, le=l + (siz << 3);
        for(; l<le; l+=4, a+=b, c+=d, e+=f){
          var m=1/a, u=c*m|0, v=e*m|0, ua=u&63, va=v&63;
          
          var j=(h << 8 & 0x3f000) + va*(h << 15 >> 19) +
                ua*((h << 8 >> 19) + va*(h >> 24));
        //   var j=0x40000;
          var i=((u ^ v)&15) + 214;
          p[l]=p[l+1]=p[l+2]=j*i >> 18;
        }
    }
};

draw= function() {
    var ms=millis(), t=ms*1e-3;
    phys();
    
    var g=new GF();
    
    var vf=g.qrot(rq, [0,0,1]),
        vr=g.qrot(rq, [1,0,0]),
        vd=g.qrot(rq, [0,1,0]);
    
    if(!loadPixels){return;}
    if(!this.imageData||!this.imageData.data||0){
        background(100);
        loadPixels();}
    if(!this.imageData||!this.imageData.data){return;}
    // var p=this.imageData.data;
    
    var vu=[100,0,0], vv=[0,100,0], vn=[0,0,1], vt=[0,0,5];
    
    var b=[];
    for(var i=0; i<4; i++){
        // b[i]=sqrt(0.25*(sin((i + 3)*ms/30) + 3));
        b[i]=(0.5*(sin((i + 3)*ms/30) + 1));
        // b[i]=0.5;
    }
    // b[0]=b[1]=b[2]=0.5;
    var h=((b[0] - b[1] - b[2] + b[3])*64 & 0xff) << 24 |
          ((b[1] - b[0])*64 & 0x7f) << 17 |
          ((b[2] - b[0])*64 & 0x7f) << 10 |
          b[0]*64 << 4;
    
    rscr(this.imageData.data,width,height,
        [vt,vn,vu,vv], [vo,vf,vr,vd],
        h);
    
    updatePixels();
    
    
    fill(250, 3, 3);
    for(var i=0, s=''; i<3; i++){s+=vo[i].toFixed(2)+"\n";}
    
    text(millis()-ms+"\n"+(this.__frameRate).toFixed(1)+"\n\n"+s+"\n"/*+vd[0].p*/, 4,180,380,380);
};
textSize(14);
