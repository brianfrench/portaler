





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
    this.abcd_adcb=function(a, b, c, d){
        return this.dot(a, b)*this.dot(c, d) -
               this.dot(a, d)*this.dot(c, b);
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

var rq=[1, 0, 0, 0], rdq=[0,0,0,0], vo=[0,0,3], vdo=[0.1, 0, 0], vhi=0;
rq=[3,-1,-1,1];
var pause=0, lighta=1;


var kp=[], keyPressed=function(){
    kp[keyCode]=1;
    // if(keyCode===SHIFT){gtex.length=0;ntex=[0,6,0,0];}
    if(keyCode==='P'.charCodeAt(0)){pause^=1;loop();}
    if(keyCode==='L'.charCodeAt(0)){lighta^=1;}
    
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
    var rm=dt*5e1;
    rq=g.norm(g.qm(rq,
        g.norm([cos(rm), sin(rm)*(kp[UP] - kp[DOWN]),
                         sin(rm)*(kp[RIGHT] - kp[LEFT]), 0])));
    
    vdo=g.mul(vdo, pow(0.1, dt));
    vdo=g.add(vdo, g.mul(g.qrot(rq, [0,0,1]),//[1,0,0],
        (kp['W'.charCodeAt(0)] - kp['S'.charCodeAt(0)])*10*dt));
    vdo=g.add(vdo, g.mul(g.qrot(rq, [1,0,0]),//[0,1,0],
        (kp['D'.charCodeAt(0)] - kp['A'.charCodeAt(0)])*10*dt));
    vdo=g.add(vdo, g.mul([0,0,-1],//g.qrot(rq, [0,1,0]), 
        (kp['Q'.charCodeAt(0)] - kp['E'.charCodeAt(0)])*10*dt));
    
    vdo=g.add(vdo, [0,0,sin(ms*2e-1)*8e-3]);
    vo=g.add(vo, g.mul(vdo, dt));
    
    if(vo[2]<0){
        vo[2]=0;
        vdo[2]=abs(vdo[2]);
    }
};

var rsc=function(face, view, tex, bri, screen){
    var g=new this.GF();
    var p=screen.p, ppl=screen.ppl;
    
    var od=g.sub(view.o, face.o);
    var zbb=g.dot(view.c, face.n),
        zbm=g.dot(view.d, face.n),
        zm= g.dot(view.r, face.n),
        
        ubb=g.abcd_adcb(od, face.u, view.c, face.n),
        ubm=g.abcd_adcb(od, face.u, view.d, face.n),
        um= g.abcd_adcb(od, face.u, view.r, face.n),
        
        vbb=g.abcd_adcb(od, face.v, view.c, face.n),
        vbm=g.abcd_adcb(od, face.v, view.d, face.n),
        vm= g.abcd_adcb(od, face.v, view.r, face.n);
    
    for(var y=0; y<screen.h; y++){
        var l=y*ppl << 2, el=l + (ppl << 2);
        var zb=zbb + y*zbm, ub=ubb + y*ubm, vb=vbb + y*vbm;
        for(; l<el; l+=4, zb+=zm, ub+=um, vb+=vm){
            var z=1/zb, u=z*ub|0, v=z*vb|0, ua=u&31, va=v&31;
            var h=bri[(u>>5&15) | (v>>1&240)];
            var b=(h << 12 & 0x7f000) + va*(h << 17 >> 17) +
                  ua*((h << 9 >> 17) + va*(h >> 21));
            // var b=(h << 12 & 0x7f000) + va*(h << 17 >> 23 << 6) +
            //       ua*((h << 9 >> 23 << 6) + va*(h >> 21));
            var t=tex[(u&63) | (v&63) << 6];
            // var t=0xffffff;
            p[l  ]=b*(t >> 16) >> 19;
            p[l+1]=b*(t >> 8 & 255) >> 19;
            p[l+2]=b*(t & 255) >> 19;
        }
    }
};



var gtex=[], ntex=[0,6,0,0], texf=function(x,y){
    if(1){
        // return[1,1,1];
    }
    var g=new GF();//this.g;
    var s=20;
    var u,v,m,n,w, m0,rn0,vx=0,vy=0;
    var h;
    for(var i=0; i<3; i++){
        var d0=1e3, d1=1e3, d2=1e3, dx,dy,d0x,d0y, rn=0x234159c;
        for(var j=0; j<18; j++){
            rn^=(rn>>3&0xfff)*(rn>>17&0xfff)+0x17e8d9cb;
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
        h=(sqrt(d1) - sqrt(d0))*24;
        if(h>1){h=1;}
        if(h<0){h=0;}
        // h*=h;
        h*=h*(3 - 2*h);
        var w=rn0*1e-8;
        m=noise(d0x, d0y, w);
        w=noise(w + m*16, d0x*6, d0y*6);
        // w*=w*(3 - 2*w);
        m=(h*0.2 + 0.8)*w;
        
        if(i===0){m0=m;x+=1e-6;}
        if(i===1){vx=m-m0;x-=1e-6;y+=1e-6;}
        if(i===2){vy=m-m0;}
    }
    var n=g.dot(g.norm(g.cross([1e-6,0,vx/10],[0,1e-6,vy/10])),
                g.norm([1,1,9]));
    // n=(n + 1)*0.5;
    
    var c0=[[0.6, 0.5, 0.4], [0.4, 0.5, 0.5]][rn0&1],
        c1=[[0.8, 0.8, 0.8], [0.8, 0.8, 0.9]][rn0&1],
        c2=[0.6, 0.6, 0.6];

    return g.mul(g.lerp(c2, g.lerp(c0, c1, m0*m0*(3 - 2*m0)), h), 0.9*(n - 1) + 1);
    // return [x,y,0.5];
};

var gbri=[];

mouseClicked=function(){
    if(mouseButton===RIGHT){
        ntex=[0,6,0,0];
    }else{
        // gbri.length=0;
    }
};

var mbri=function(lp, ld){
    var g=new this.GF();
    var br=[];
    // var lp=[0.1,0.1,1], ld=g.norm([1,1,-3]);
    // lp=[sin(t*50)*0.6 + 0.5,
    //     cos(t*40)*0.6 + 0.5,
    //     0.3 + sin(t*20)*0.2];
    
    // lp=[0.5,0.5,0.2];
    
    for(var i=0; i<256; i++){
        var x=(i&15)/16, y=(i>>4)/16;
        
        var n=0;
        var d=g.sub([x,y,0], lp);
        var m=g.dot(g.norm(d), ld);
        // m=(m*5-4.6);
        // m=(m*20 - 19);
        // m=m*0.2 - 0.13;
        // m=pow
        // m-=0.9;
        m=m*10 - 9;
        if(m<0){m=0;}
        m=m*0.5 + 1*this.pow(1-4*this.pow(m-0.5,2),10) +
            2*this.pow(m,16);
        // m*=m;
        m*=0.2/g.mag2(d);
        m=m>1?1:m<0?0:m;
        n+=m;
        // }
        n=n>1?1:n<0?0:n;
        br[i]=0.95*(n - 1) + 1;
        
    }
    for(var j=0; j<1; j++){
        var cr=br.slice();
        for(var i=0; i<256; i++){
            var x=i&15, y=i&240;
            // br[i]=(cr[i]*2 +
            //         cr[(x+1&15)|y] +
            //         cr[x|(y+16&240)] +
            //         cr[(x-1&15)|y] +
            //         cr[x|(y-16&240)])/6;
            var m=0;
            for(var k=0; k<9; k++){
                m+=cr[(x+k%3&15)|(y+(k/3|0)*16&240)]*
                    // [53,97,53, 97,100,97, 53,97,53][k];
                    [9,46,9, 46,100,46, 9,46,9][k];
            }
            br[i]=m/320;
        }
    }
    for(var i=0; i<256; i++){br[i]=this.pow(br[i], 0.45)*127.9|0;}
    var bri=[];
    for(var i=0; i<256; i++){
        // var 
        var x=i&15, y=i&240,
            a=br[i], b=br[(x+1&15)|y],
            c=br[x|(y+16&240)], d=br[(x+1&15)|(y+16&240)];
        var h=(a - b - c + d & 0x1ff) << 23 |
              (b - a & 0xff) << 15 |
              (c - a & 0xff) << 7 |
              a;
        bri[i]=h;
    }
    //gbri=bri.slice();
    // cos(0);
    return bri;
};

draw= function() {
    if(pause){noLoop();return;}
    if(gtex.length<0x4000){
        for(var i=0; i<0x4000; i++){
            gtex[i]=(random(16777216)&0x3f3f3f)+0xa0a0a0;
        }
    }
    
    var g=new GF();
    
    var ms=millis(), t=ms*1e-3;
    
    // if(gtex.length!==0x1000){
    //     gtex=[];
    //     for(var i=0; i<0x1000; i++){
    //         gtex[i]=((i&63)|(i>>6))*0x20202 + 0x305030;
    //     }
    //     // return;
    // }
    
    while(ntex[1]>=0){
        if(((ntex[0]&7)===0||ntex[1]>4)&&millis()-ms>6){break;}
        var s=1<<ntex[1];
        if(((ntex[2]|ntex[3])&s) || ntex[1]===6){
            var vc=[0,0,0], x=ntex[2]/64, y=ntex[3]/64;
            for(var i=0; i<4; i++){
                vc=g.add(vc, texf(x+(i&1)/128, y+(i>>1)/128));
            }
            vc=g.mul(vc, 0.25);
            var c=color(sqrt(vc[0])*256, sqrt(vc[1])*256, sqrt(vc[2])*256)&0xffffff;
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
    
    
    
    var vf=g.qrot(rq, [0,0,1]),
        vr=g.qrot(rq, [1,0,0]),
        vd=g.qrot(rq, [0,1,0]);
    
    if(!loadPixels){return;}
    if(!this.imageData || !this.imageData.data || 0){
        background(0xff552255); loadPixels();}
    if(!this.imageData || !this.imageData.data){return;}
    
    var vu=[100,0,0], vv=[0,100,0], vn=[0,0,1], vt=[0,0,0];
    
    if(!gbri.length||lighta){
        var m=mouseX*0.5 + t*30;
        var lp=[0.5*(sin(m)*0.7 + sin(m*3)*0.5 + 1),
                0.5*(cos(m)*0.7 - cos(m*3)*0.5 + 1),
                0.4];
                // 0.1 + (1 - 4*sq(t*0.5%1 - 0.5))*0.2];
        // var lp=[0.5*(sin(t*50)*0.7 + sin(t*150)*0.5 + 1),
        //         0.5*(cos(t*50)*0.7 - cos(t*150)*0.5 + 1),
        //         0.3];
        // var ld=g.norm(g.sub([mouseX/width,mouseY/height,0.0], lp));
        var ld=g.norm(g.sub([0.5,0.5,lp[2]-0.3], lp));
        ld=g.qrot([cos(sin(m*4)*10), 0, 0, sin(sin(m*4)*10)], ld);
        // var ld=g.norm(g.sub([0.5,0.5,0], lp));
        // var ld=[0,0,-1];
        
        gbri=mbri(lp, ld);
    }
    
    // rscr(this.imageData.data, width, height,
    //     [vt,vn,vu,vv], [vo,vf,vr,vd],
    //     gtex, bri);
        
    rsc({p:[[]], o:vt, n:vn, u:vu, v:vv},
        {   o:vo,//[1,1,-3],
            c:g.sub(vf, g.add(g.mul(vr, width/height), vd)),
            r:g.mul(vr, 2/height),
            d:g.mul(vd, 2/height)},
        gtex,
        gbri,
        {p:this.imageData.data,
            x:0, y:0, w:width, h:height, ppl:width});
    
    updatePixels();
    
    
    fill(250, 3, 3);
    for(var i=0, s=''; i<3; i++){s+=vo[i].toFixed(2)+"\n";}
    
    text(millis()-ms+"\n"+(this.__frameRate).toFixed(1)+"\n\n"+s+"\n"/*+vd[0].p*/, 4,16);
};
textSize(14);

