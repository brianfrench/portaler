
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
frameRate(30);
var rsc=function(face, view, tex, screen){
    var g=new this.GF();
    var p=screen.p, ppl=screen.ppl;
    
    var od=g.sub(view.o, face.o);
    var zbb=g.dot(view.c, face.n),
        zbm=g.dot(view.d, face.n),
        zm=g.dot(view.r, face.n),
        
        ubb=g.abcd_adcb(od, face.u, view.c, face.n),
        ubm=g.abcd_adcb(od, face.u, view.d, face.n),
        um=g.abcd_adcb(od, face.u, view.r, face.n),
        
        vbb=g.abcd_adcb(od, face.v, view.c, face.n),
        vbm=g.abcd_adcb(od, face.v, view.d, face.n),
        vm=g.abcd_adcb(od, face.v, view.r, face.n);
    
    for(var y=0; y<screen.h; y++){
        var l=y*ppl << 2, el=l + (ppl << 2);
        var zb=zbb + y*zbm, ub=ubb + y*ubm, vb=vbb + y*vbm;
        for(; l<el; l+=4, zb+=zm, ub+=um, vb+=vm){
            var m=1/zb, u=m*ub, v=m*vb;
            p[l]=p[l+1]=p[l+2]=(u^v)&255;
        }
    }
    
    // for(var y=0; y<50; y++){
    //     var l=y*ppl << 2;
    //     for(var x=0; x<50; x++, l+=4){
    //         // p[l]=p[l+1]=p[l+2]=view.c[0]*255&255;
    //         var vv=g.add(g.add(g.mul(view.r, x), g.mul(view.d, y)),
    //             view.c);
    //         var a=g.dot(g.sub(face.o, view.o), face.n)/
    //             g.dot(vv, face.n);
    //         var vp=g.add(view.o, g.mul(vv, a));
    //         p[l]=p[l+1]=p[l+2]=((vp[0]*20 ^ vp[1]*20)&63) + 128;
            
    //     }
    // }
    
    // for(var l=screen.h*screen.ppl<<2; l>=0; l-=4){
    //     p[l]=p[l+1]=p[l+2]=l&255;
    // }
};

var viewer={o:[1,1,-3], d:[0,0,0]};

draw= function() {
    var ms=-millis();
    
    var g=new this.GF();
    
    if(!loadPixels){return;}
    if(!this.imageData || !this.imageData.data || 0){
        background(0xff552255); loadPixels();}
    if(!this.imageData || !this.imageData.data){return;}
    
    var q=[1,0,0,0], i;
    for(i=0; i<7; i++){
        q=g.qm([1, (height/2 - mouseY)*1e-3, (mouseX - width/2)*1e-3, 0], q);
    }
    q=g.norm(q);
    
    viewer.o=g.add(viewer.o, viewer.d);
    viewer.d=g.add(g.mul(viewer.d, 0.7), 
        g.qrot(q, [0,0,5e-2*(mouseIsPressed*(mouseButton===LEFT?1:-1))]));
    
    
    rsc({p:[[]], o:[0,0,0], n:[0,0,1], u:[64,0,0], v:[0,64,0]},
        {   o:viewer.o,//[1,1,-3],
            c:g.qrot(q, [-1,-1,1]),
            r:g.qrot(q, [2/width,0,0]),
            d:g.qrot(q, [0,2/height,0])},
        {},
        {p:this.imageData.data, x:0, y:0, w:width, h:height, ppl:width}         );
    
    updatePixels();
    
    ms+=millis();
    fill(255, 0, 0);
    text(ms+"\n"+this.__frameRate.toFixed(1), 4, height-24);
};
