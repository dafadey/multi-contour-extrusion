defaultpen(fontsize(8pt));
unitsize(22cm);

struct interval {
  real[] reds; // to be removed
  real[] greens; // to be added
  real[] added;
  int[] left={0,0};
  int[] right={0,0};
};

/*
interval a0;
interval a1;
interval a2;
interval a3;
interval a4;
interval a5;
interval a6;
interval a7;
a7.reds = new real[]{0.5000013038471354, };
interval a8;
interval a9;
interval a10;
interval a11;
interval a12;
a12.reds = new real[]{0.5, };
interval a13;
interval a14;
a14.reds = new real[]{0.4999967445317358, };
interval a15;
interval a16;
a16.reds = new real[]{0.5000000926662145, };
interval a17;
interval a18;
interval a19;
interval a20;
interval a21;
interval a22;
interval a23;
interval a24;
interval a25;
interval a26;
interval a27;
interval a28;
interval a29;
interval a30;
interval a31;
interval a32;
a32.reds = new real[]{0.4999996423709795, };
interval a33;
a33.reds = new real[]{0.4999985173406054, };
interval a34;
interval a35;
interval a36;
interval a37;
interval a38;
interval a39;
interval a40;
interval a41;
interval a42;
interval a43;
a43.reds = new real[]{0.49999987240939464, };
interval a44;
interval a45;
interval a46;
interval a47;
interval a48;
interval a49;
interval a50;
interval a51;
interval a52;
interval a53;
interval a54;
interval a55;
interval a56;

interval b0;
b0.added = new real[]{0.570595174521068, };
interval b1;
interval b2;
b2.added = new real[]{0.45065591548828043, };
interval b3;
interval b4;
interval b5;
interval b6;
interval b7;
interval b8;
interval b9;
interval b10;
interval b11;
interval b12;
interval b13;
interval b14;
interval b15;
interval b16;
interval b17;
interval b18;
interval b19;
interval b20;
interval b21;
interval b22;
interval b23;
b23.added = new real[]{0.46656478000514, };
interval b24;
interval b25;
interval b26;
interval b27;
interval b28;
b28.added = new real[]{0.5031523701817604, };
interval b29;
interval b30;
interval b31;
interval b32;
interval b33;
interval b34;
interval b35;
interval b36;
interval b37;
interval b38;
interval b39;
interval b40;
interval b41;
interval b42;
interval b43;
interval b44;
interval b45;
interval b46;
interval b47;
interval b48;
interval b49;
interval b50;
interval b51;
b51.added = new real[]{0.5028323194205925, };
interval b52;
interval b53;
interval b54;
interval b55;
interval b56;
*/



interval a0;
a0.added=new real[]{0.11505818016928723, 0.2305120460375735, 0.3683522074322502, 0.5179902743309245, 0.5937867584100248, 0.7404925309401128, 0.911570965167776};

interval a1;
a1.added=new real[]{0.09853226867691238, 0.2854253535525984, 0.3705512565647206, 0.5219047142530894, 0.6653851339751367, 0.7491859475521024, 0.8787305964649722};

interval a2;
a2.added=new real[]{0.11083978888103067, 0.22457691797506046, 0.3434197497357221, 0.4678841460922509, 0.5973259254985811, 0.7303767583097722, 0.8652562817075624};

interval a3;
a3.added=new real[]{0.1251482279939197, 0.24665335204470615, 0.3854537828560535, 0.5301602699583904, 0.6037391179612567, 0.7468081240429763, 0.873003638617264};

interval b0;
b0.added=new real[]{0.06857284744251976, 0.11505818016928723, 0.1693704667360144, 0.2305120460375735, 0.29727544258935457, 0.3683522074322502, 0.44238384618476334, 0.5179902743309245, 0.5937867584100248, 0.6684013716249498, 0.7404925309401128, 0.824437869308812, 0.911570965167776};
//b0.added=new real[]{0.5179902743309245};

interval b1;
b1.added=new real[]{0.09853226867691238, 0.19428804507883005, 0.2854253535525984, 0.3705512565647206, 0.44909159305573804, 0.5219047142530894, 0.5920897507803214, 0.6653851339751367, 0.7491859475521024, 0.8787305964649722};
//b1.added=new real[]{0.5219047142530894};

interval b2;
b2.added=new real[]{0.11083978888103067, 0.22457691797506046, 0.3434197497357221, 0.4678841460922509, 0.5973259254985811, 0.7303767583097722, 0.8652562817075624};
//b2.added=new real[]{0.4678841460922509};

interval b3;
b3.added=new real[]{0.1251482279939197, 0.24665335204470615, 0.3155244083943529, 0.3854537828560535, 0.4570920053087017, 0.5301602699583904, 0.6037391179612567, 0.6764901560856981, 0.7468081240429763, 0.8129271113247586, 0.873003638617264, 0.9252078343525193};
//b3.added=new real[]{0.5301602699583904};


real[] get_params_for_interval(interval i) {
  real[] res={.0};
  for(real v : i.reds)
    res.push(v);
  for(real v : i.added)
    res.push(v);
  res.push(1.);
  return res;
}

int[][] connect_intervals(interval[] intervals0, interval[] intervals1) {
  int[][] indexes;
  if(intervals0.length !=intervals1.length)
    return indexes;
  int i0_stride=0;
  for(interval i : intervals0) 
    i0_stride += 1+i.reds.length + i.added.length;
  int i1_stride=0;
  for(interval i : intervals1) 
    i1_stride += 1+i.reds.length + i.added.length;
  write("i0_stride=", i0_stride);
  write("i1_stride=", i1_stride);
  int i0=0;
  int i1=0;
  for(int i=0;i<intervals0.length;++i) {
    int ii1=0;
    int ii0=0;
    write("i0="+string(i0)+", i1="+string(i1));
    real[] inter0ext = get_params_for_interval(intervals0[i]);
    real[] inter1ext = get_params_for_interval(intervals1[i]);
    while(ii0 < inter0ext.length - 1 || ii1 < inter1ext.length - 1) {
      //write("ii0="+string(ii0)+", ii1="+string(ii1)+" 0ext_len="+string(inter0ext.length)+" 1ext_len="+string(inter1ext.length));
      real dx11 = (ii0+1 == inter0ext.length || ii1+1 == inter1ext.length) ? 1.0 : abs(inter0ext[ii0+1] - inter1ext[ii1+1]);
      real dx10 = ii0+1 == inter0ext.length ? 1.0 : abs(inter0ext[ii0+1] - inter1ext[ii1]);
      real dx01 = ii1+1 == inter1ext.length ? 1.0 : abs(inter0ext[ii0] - inter1ext[ii1+1]);
      write("dx11="+string(dx11)+" dx10="+string(dx10)+" dx01="+string(dx01));
      if(dx11 <= dx10 && dx11 <= dx01) {
        int[] q;
        q.push(i0+ii0);
        q.push((i0+ii0+1) % i0_stride);
        q.push(i0_stride+(i1+ii1+1) % i1_stride);
        q.push(i0_stride+i1+ii1);
        indexes.push(q);
        ii0+=1;
        ii1+=1;
      }
      if(dx10 <= dx11 && dx10 <= dx01) {
        int[] t;
        t.push(i0+ii0);
        t.push((i0+ii0+1) % i0_stride);
        t.push(i0_stride+(i1+ii1) % i1_stride);
        indexes.push(t);
        ii0+=1;
      }
      if(dx01 <= dx11 && dx01 <= dx10) {
        int[] t;
        t.push((i0+ii0) % i0_stride);
        t.push(i0_stride+(i1+ii1+1) % i1_stride);
        t.push(i0_stride+i1+ii1);
        indexes.push(t);
        ii1+=1;
      }
    }
    i0 += inter0ext.length-1;
    i1 += inter1ext.length-1;
  }
  return indexes;
}

interval[] ias = {a0, a1, a2, a3};
interval[] ibs = {b0, b1, b2, b3};

/*
interval[] ias = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20,  a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56};
interval[] ibs = {b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20,  b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b36, b37, b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48, b49, b50, b51, b52, b53, b54, b55, b56};
*/

real[] make_pts(interval[] ias) {
  real[] a;

  for(int i=0; i<ias.length; ++i) {
    interval inter = ias[i];
    int n = ias.length;
    real x0=i/n;
    real x1=(i+1)/n;
    a.push(x0);
    for(real r : inter.reds)
      a.push(x0+r*(x1-x0));
    for(real r : inter.added)
      a.push(x0+r*(x1-x0));
  }
  return a;
}

real[] a = make_pts(ias);
real[] b = make_pts(ibs);

//write(a);
write(a.length);
//write(b);  
write(b.length);
int[][] indexes = connect_intervals(ias, ibs);

write(indexes);

int stride = a.length;
real x=0;
for(int id=0; id<indexes.length; ++id) { 
  int[] i = indexes[id];
  guide g;
  for(int j=0;j<i.length;++j)
  {
    if(i[j]>=stride+b.length)
      continue;
    pair p=(i[j]<stride ? a[i[j]] : b[i[j]-stride], i[j]<stride ? 0. : .3);
    g=g--p;
  }
  x+=0.003;
  if(id == indexes.length - 1)
    draw(shift(x)*(g--cycle));
  else
    filldraw(shift(x)*(g--cycle),red);
}
