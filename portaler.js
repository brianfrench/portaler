




// click, right-click

frameRate(30);
mouseX=width>>1|7;
mouseY=height>>1|7;

var world=[
    [
        {p:[[0,0,3],[0,3,3],[0,3,0],[0,0,0]], i:-1},
        {p:[[4,0,3],[4,3,3],[0,3,3],[0,0,3]], i:-1},
        {p:[[0,3,3],[4,3,3],[4,3,0],[0,3,0]], i:-1},
        {p:[[0,3,0],[4,3,0],[4,0,0],[0,0,0]], i:-1},
        {p:[[0,0,0],[4,0,0],[4,0,3],[0,0,3]], i:-1},
        {p:[[4,3,0],[4,3,3],[5,2,2],[5,2,1]], i:-1},
        {p:[[4,3,3],[4,0,3],[5,1,2],[5,2,2]], i:-1},
        {p:[[4,0,3],[4,0,0],[5,1,1],[5,1,2]], i:-1},
        {p:[[4,0,0],[4,3,0],[5,2,1],[5,1,1]], i:-1},
        {p:[[5,2,1],[5,2,2],[5,1,2],[5,1,1]], i:1}
    ],
    [
        // {p:[[6,0,0],[6,3,0],[9,3,3],[9,0,3]], i:-1},
        {p:[[5,2,1],[5,2,2],[5.5,2,2],[5.5,2,1]], i:-1},
        {p:[[5,1,2],[5,1,1],[5.5,1,1],[5.5,1,2]], i:-1},
        {p:[[5,1,1],[5,2,1],[5.5,2,1],[5.5,1,1]], i:-1},
        {p:[[5,2,2],[5,1,2],[5.5,1,2],[5.5,2,2]], i:-1},
        {p:[[5,1,2],[5,2,2],[5,2,1],[5,1,1]], i:0},
        {p:[[5.5,1,1],[5.5,2,1],[5.5,2,2],[5.5,1,2]], i:2}
    ],
    [
        {p:[[6,0,0],[6,3,0],[9,3,3],[9,0,3]], i:-1},
        {p:[[5.5,1,2],[5.5,2,2],[5.5,2,1],[5.5,1,1]], i:1}
        
    ]
];

var textures=[
    
];
// a

var gf=function(){
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
                // var m=m1/(m1 - m0);
                op[j++]=this.add(v0,
                    this.mul(this.sub(v1, v0), m0/(m0 - m1)));
            }
            if(m1>=0){
                op[j++]=v1.slice();
            }
            m0=m1;
            v0=v1;//.slice();
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

var rpoly=function(p3, pt, po, p,wid,hei,tex,br){
    var g=new this.gf();//this.g;
    var i;
    var p2=[], v, n=p3.length;
    if(n<=0){return;}
    for(i=0; i<n; i++){
        v=g.sub(p3[i], po[0]);
        var m=0.49/g.dot(v, po[1]);
        p2[i]=[wid*(0.5 + g.dot(v, po[2])*m),
               hei*(0.5 + g.dot(v, po[3])*m)];
    }
    
    for(var y=hei, i=0, j=0; i<n; i++){
        if(p2[i][1]<y){
            j=i;
            y=p2[i][1];
        }
    }
    
    y=y + 1|0;
    if(y<0){y=0;}
    var i0=j - 1, pi0=j, i1=j + 1, pi1=j;
    
    var vp=g.sub(po[0], pt[0]),
        v0=g.add(g.add(g.mul(po[2], -1/0.98),
            g.mul(po[3], -1/0.98)), po[1]);
    var num=-g.dot(vp, pt[1]),
        ou=g.dot(vp, pt[2]), ov=g.dot(vp, pt[3]),
        a0=g.dot(v0, pt[1]),
            du_da0=g.dot(v0, pt[2]), dv_da0=g.dot(v0, pt[3]),
        da_dx=g.dot(po[2], pt[1])/(wid*0.49),
        da_dy=g.dot(po[3], pt[1])/(hei*0.49),
        du_dx=g.dot(po[2], pt[2])/(wid*0.49),
        du_dy=g.dot(po[3], pt[2])/(hei*0.49),
        dv_dx=g.dot(po[2], pt[3])/(wid*0.49),
        dv_dy=g.dot(po[3], pt[3])/(hei*0.49);
    
    // var b0=br[0]*0x100000|0,
    //     b1=(br[1] - br[0])*0x4000|0,
    //     b2=(br[0] - br[1] - br[2] + br[3])*0x100|0,
    //     b3=(br[2] - br[0])*0x4000|0;
    
    var b=br*256|0;
    
    while(y<hei){
        if(i0<0){i0=n-1;}
        if(i1>=n){i1=0;}
        var y0=p2[pi0][1], y1=p2[i0][1],
            y2=p2[pi1][1], y3=p2[i1][1],
            x0=p2[pi0][0], x1=p2[i0][0],
            x2=p2[pi1][0], x3=p2[i1][0];
        if(y1<y0 || y3<y2){break;}
        if(y1<y){pi0=i0; i0--; continue;}
        if(y3<y){pi1=i1; i1++; continue;}
        
        var xl=(x1 - x0)*(y - y0)/(y1 - y0) + x0 + 1|0,
            xr=(x3 - x2)*(y - y2)/(y3 - y2) + x2 + 1|0;
        if(xl>xr){var t=xl; xl=xr; xr=t;}
        if(xl<0){xl=0;}
        if(xr>wid){xr=wid;}
        // xl=xl+1|0;xr=xr+1|0;
        var l=y*wid + xl << 2;
        
        
        var a=a0 + da_dx*xl + da_dy*y,
            du=du_da0 + du_dx*xl + du_dy*y,
            dv=dv_da0 + dv_dx*xl + dv_dy*y;
        
        // var mm=1/a;
            // mm*=2 - mm*a;
            // var m=num*mm;
        for(var el=l + (xr - xl << 2); l<el;
            l+=4, du+=du_dx, dv+=dv_dx, a+=da_dx){
            var m=num/a;
            var c=tex[((ou + m*du)&0x3f)|((ov + m*dv)&0x3f)<<6];
            //     b=b0 + u*(b1 + v*b2) + v*b3;
            
            p[l  ]=(c>>16)*b >> 8;
            p[l+1]=(c>>8&0xff)*b >> 8;
            p[l+2]=(c&0xff)*b >> 8;
        }
        
        y++;
    }
    // sin(0);
};

// var smax=function(x,y,m){
//     m=(x - y)*m;
//     if(m>1){return 1;}
//     if(m<-1){return 0;}
//     return 0.5 + m*(0.75 - 0.25*m*m);
// };

var gtex=[], ntex=[0,6,0,0], texf=function(x,y){
    var g=new gf();
    var s=20;
    var u,v,m,n,w, m0,rn0,vx=0,vy=0;
    var h;
    for(var i=0; i<3; i++){
        var d0=1e3, d1=1e3, d2=1e3, dx,dy,d0x,d0y, rn=0x314159c;
        for(var j=0; j<9; j++){
            rn^=(rn&0xfff)*(rn>>17&0xfff)+0x17e8d9cb;
            rn^=rn<<j+9;
            for(var k=0; k<9; k++){
                dx=x - (rn&1023)/1024 + (k%3) - 1;
                dy=y - (rn>>10&1023)/1024 + ~~(k/3) - 1;
                m=dx*dx + dy*dy;
                if(m<d0){
                    d2=d1; d1=d0; d0=m; d0x=dx; d0y=dy; rn0=rn;
                }else if(m<d1){
                    d2=d1; d1=m;
                }else if(m<d2){
                    d2=m;
                }
            }
        }
        h=(sqrt(d1) - sqrt(d0))*15;
        if(h>1){h=1;}
        if(h<0){h=0;}
        h*=h*(3 - 2*h);
        var w=rn0*1e-8;
        m=noise(d0x, d0y, w);
        w=noise(w + m*16, d0x*6, d0y*6);
        w*=w*(3 - 2*w);
        m=(h*0.5 + 0.5)*w;
        
        if(i===0){m0=m;x+=1e-6;}
        if(i===1){vx=m-m0;x-=1e-6;y+=1e-6;}
        if(i===2){vy=m-m0;}
    }
    var n=g.dot(g.norm(g.cross([1e-6,0,vx/10],[0,1e-6,vy/10])),
                g.norm([1,1,9]));
    n=(n + 1)*0.5;
    // var 
    var c0=[0xcc9955, 0x778899][rn0&1],
        c1=[0xeeee99, 0xeeeeee][rn0&1],
        c2=180;//lerpColor(c0,c1,0.5);
    return lerpColor(0, lerpColor(c2, lerpColor(c0, c1, m0*m0*(3 - 2*m0)), h), 0.5*(n - 1) + 1)&0xffffff;
};

var kp=[], keyPressed=function(){
    kp[keyCode]=1;
    if(keyCode===SHIFT){gtex.length=0;ntex=[0,6,0,0];}
    if(keyCode==='M'.charCodeAt(0)){mapon^=1;}
    
},keyReleased=function(){kp[keyCode]=0;};
for(var i=0; i<256; i++){kp[i]=0;}

var view=function(world, hi, o, beam){
    var g=new this.gf();
    var vd=[], stack=[{hi:hi, o:o, beam:beam, g:1}];
    while(stack.length){
        var so=stack[stack.length - 1];
        if(!so.g){stack.length--; continue;}
        var h=world[so.hi];
        for(var i=0; i<h.length; i++){
            var p=h[i];
            if(g.dot(g.sub(so.o, p.p[0]),
                g.cross(g.sub(p.p[0], p.p[1]),
                        g.sub(p.p[1], p.p[2]))) > 0){
                continue;
            }
            var cp=g.cull(p.p, so.o, so.beam);
            if(p.i>=0){
                stack[stack.length]={hi:p.i, o:o,
                    beam:g.beamify(cp, so.o), g:1};
            }else{
                vd[vd.length]={
                    // p:[[0,0,0],[1,1,1],[2,2,2]],//g.cull(p, so.o, so.beam),
                    p:cp, hi:so.hi, pi:i
                };
            }
        }
        so.g=0;
        // stack.length--;
    }
    return vd;
};

var lms=0;
var phys=function(){
    var g=new gf();
    
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
    
    for(var i=-1; i<world.length; i++){
        var h=world[i>=0?i:vhi];
        for(var j=0; j<h.length; j++){
            var p=h[j].p;
            if(g.dot(g.sub(p[0], vo),
                g.cross(g.sub(p[1], p[0]), g.sub(p[2], p[1])))<0){
                j=0;
                break;
            }
        }
        if(j){
            vhi=i>=0?i:vhi;
            break;
        }
    }
};

draw= function() {
    var ms=millis(), t=ms*1e-3;
    
    if(gtex.length!==0x1000){
        gtex=[];
        for(var i=0; i<0x1000; i++){
            gtex[i]=((i&63)|(i>>6))*0x20202 + 0x305030;
        }
        // return;
    }
    
    while(ntex[1]>=0){
        if(((ntex[0]&7)===0||ntex[1]>4)&&millis()-ms>3){break;}
        // gtex[ntex]=texf((ntex&127)/128, (ntex>>7)/128);
        var s=1<<ntex[1];
        if(((ntex[2]|ntex[3])&s) || ntex[1]===6){
            var c=texf(ntex[2]/64, ntex[3]/64);
            for(var y=ntex[3]; y<ntex[3]+s; y++){
                for(var l=y*64 + ntex[2], i=0; i<s; i++, l++){
                    gtex[l]=c;
                }
            }
            ntex[0]++;
        }
        ntex[2]+=s;
        if(ntex[2]>=64){
            ntex[2]=0;
            ntex[3]+=s;
            if(ntex[3]>=64){
                ntex[3]=0;
                ntex[1]--;
            }
        }
        if(ntex[0]<10){break;}
    }
    
    phys();
    
    var g=new gf();
    
    var vf=g.qrot(rq, [0,0,1]),
        vr=g.qrot(rq, [1,0,0]),
        vd=g.qrot(rq, [0,1,0]);
    
    var beam=[g.sub(vf, vr), g.add(vf, vr),
              g.sub(vf, vd), g.add(vf, vd)];
    // var beam=g.beamify([[0.5,0.5,0.5],[1,0.5,0.5],[1,1,1]], vo);
    // var beam=g.beamify(
    //     g.cull([[0.5,0.5,0.5],[1,0.5,0.5],[1,1,1]], vo, beam), vo);
    
    
    if(!loadPixels){return;}
    if(!this.imageData||!this.imageData.data||0){
        background(100);
        loadPixels();}
    if(!this.imageData||!this.imageData.data){return;}
    // var p=this.imageData.data;
    
    // background(100);
    var vdat=view(world, vhi, vo, beam);
    
    for(var i=0; i<vdat.length; i++){
        var face=world[vdat[i].hi][vdat[i].pi].p,
            poly=vdat[i].p;
        
        var v0=g.sub(face[1], face[0]), v1=g.sub(face[2], face[1]),
            vn=g.norm(g.cross(v0, v1));
        var vu=g.cross(v1, vn), vv=g.cross(/*v0*/vu, vn);
        // vu=g.mul(vu, 2/g.dot(vu, v0));
        // vv=g.mul(vv, 2/g.dot(vv, v1));
        vu=g.mul(g.norm(vu), 32);
        vv=g.mul(g.norm(vv), 32);
        
        var b=sqrt(g.dot(vn, g.norm([5, 9, 3]))*0.4 + 0.6);
        
        rpoly(vdat[i].p,
            [face[0], vn, vu, vv],
            [vo, vf, vr, vd],
            this.imageData.data, width, height, gtex, b);
            
        // rpoly(vdat[i].p,
        //     face[0],vn,vu,vv,
        //     vo,vf,vr,vd,
        //     this.imageData.data,width,height,gtex,b);
        
        // noFill();stroke(255);fill(0x33777777);
        // beginShape();
        // // vertex(0,0);
        // for(var k=0; k<poly.length; k++){
        //     var v=g.sub(poly[k], vo), m=1/g.dot(v, vf);
        //     vertex((m*g.dot(v, vr) + 1.02)*(width*0.49),
        //           (m*g.dot(v, vd) + 1.02)*(height*0.49));
        //     // vertex((v[0] + 1.02)*(400*0.49),
        //     //       (v[1] + 1.02)*(400*0.49));
        // }
        // // vertex(400,400);
        // endShape(CLOSE);
    }
    
    updatePixels();
    
    fill(0x33ffffff);noFill();
    for(var i=0; mapon && i<world.length; i++){
        stroke(i===vhi?0xbbff9999:0x77ffffff);
        var h=world[i];
        for(var j=0; j<h.length; j++){
            // var p=h[j].p;
            var p=g.pclip(h[j].p, g.add(vo, g.mul(vf, -4)), vf);
            beginShape();
            for(var k=0; k<p.length; k++){
                var v=g.sub(p[k], vo), m=1/(g.dot(v, vf) + 8);
                vertex((m*g.dot(v, vr) + 1.02)*(width*0.49),
                        (m*g.dot(v, vd) + 1.02)*(height*0.49));
            }
            endShape(CLOSE);
        }
    }
    
    fill(255);
    for(var i=0, s=''; i<3; i++){s+=vo[i].toFixed(2)+"\n";}
    // s+="\n";
    // for(var i=0; i<vdat.length; i++){
    //     var p=vdat[i].p;
    //     for(var j=0; p && j<p.length; j++){
    //         var v=p[j];
    //         s+='[';
    //         for(var k=0; k<v.length; k++){
    //             s+=v[k].toFixed(1)+" ";
    //         }
    //         s+='],';
    //     }
    //     s+="\n\n";
    // }
    
    text(millis()-ms+"\n"+(this.__frameRate).toFixed(1)+"\n\n"+s+"\n"/*+vd[0].p*/, 4,8,380,380);
    
    if(ntex[0]<0x1000){
        noStroke();
        fill(0x77ffffff);
        arc(30,60,20,20,0,ntex[0]/0x1000*360);
    }
};
textSize(14);
