unitsize(22cm);
defaultpen(fontsize(8pt));
int n1,n2;

/*
n1 = 17;
n2 = 17;
srand(11);
*/

//error in last interval
n1 = 7;
n2 = 17;
srand(11);


//real[] reds={0.0, 0.25, 0.5, 0.75, 1.0};
//real[] greens={0.0, 0.007856472074541961, 0.017699564031766573, 0.02969804672311275, 0.04371677030908083, 0.05949822518826462, 0.07673074587711311, 0.09507660429523074, 0.11418515497303483, 0.1337001797400589, 0.15326426047136238, 0.1725232846122233, 0.1911309717363811, 0.21279838015004732, 0.23528859117321552, 0.25811330128303317, 0.28077624339021157, 0.30280057761254, 0.323762639823019, 0.3433420470333703, 0.3614067391679919, 0.3781541415418215, 0.39429707056360663, 0.41115539574349536, 0.43043002538646663, 0.4602259787125168, 0.4881185819032354, 0.5156616343231812, 0.5439246591841672, 0.5734564219019046, 0.6043851112098962, 0.6365506516851341, 0.6696130212514412, 0.7031298099500385, 0.7366128518581251, 0.7695752867244527, 0.8015781748328141, 0.819717925954266, 0.8381364244650785, 0.8570050115748348, 0.876250253420592, 0.8956299763596315, 0.9147916648155889, 0.9333125138589882, 0.9507274053571657, 0.9665507905510062, 0.9803007047867554, 0.9915409675743779, 1.0};


/*
real[] reds={0.0, 0.015624978468736032, 0.031249968113359533, 0.04687503924858752, 0.06250002889321103, 0.07812500736194705, 0.09374998583068309, 0.10937497547530659, 0.12500004661053457, 0.1406250362551581, 0.1562500147238941, 0.17187500436851763, 0.1875000866795668, 0.20312505397240738, 0.21875004361703088, 0.23437503326165438, 0.25000000055449495, 0.26562482995242587, 0.2812498628162265, 0.2968749115415806, 0.31249995019703286, 0.3281249830608335, 0.3437500317861876, 0.3593748453223424, 0.3749998940476965, 0.3906249327031488, 0.40624997563685494, 0.4218750142923072, 0.43750004715610785, 0.4531248765540388, 0.46874991520949105, 0.48437495814319725, 0.49999999679864954, 0.515624986443273, 0.5312499537361136, 0.5468750360471628, 0.5625000256917864, 0.5781249929846269, 0.5937499826292505, 0.609374972273874, 0.6250000545849231, 0.6406250218777638, 0.6562500115223873, 0.6718750011670108, 0.6875000723022389, 0.7031250507709749, 0.7187500404155984, 0.7343750188843345, 0.750000008528958, 0.765624830010555, 0.7812498657556297, 0.7968749094751375, 0.8124999502842687, 0.8281249910934, 0.8437500348129077, 0.8593748562945046, 0.8749998971036359, 0.8906249328487107, 0.9062499736578419, 0.9218750173773496, 0.9375000581864809, 0.9531248796680778, 0.9687499208846617, 0.984374961693793, 1.0};

real[] greens={0.0, 0.0112403485019614, 0.01969934027464576, 0.027555871925871483, 0.037398905616329095, 0.04939744407031521, 0.06341614911991533, 0.07919760993121953, 0.09643013166402281, 0.11477594646242882, 0.13388456836221635, 0.15339950199242097, 0.1729636312556734, 0.192222680105898, 0.21083033890488953, 0.23249773894487727, 0.25498799008382905, 0.2778126949487838, 0.300475577786921, 0.3224998974269634, 0.3434620353038058, 0.36304143573421177, 0.38110611147022566, 0.39785354621670005, 0.41399639946563505, 0.4224256469101172, 0.4308548629239534, 0.4404921970979461, 0.4501295276926906, 0.46502739288819667, 0.4799252580837028, 0.4938715989842561, 0.5078179398848094, 0.5215895264492808, 0.5353608816701554, 0.5494924470734346, 0.5636239716585927, 0.5931557505588854, 0.6240844286874968, 0.640167265156593, 0.6562500152378353, 0.6727812500801675, 0.6893124260333412, 0.7060707983583356, 0.72282917068333, 0.7395707218345998, 0.7563122293112702, 0.7727934504571123, 0.7892746831414387, 0.8052761193775502, 0.8212775259127844, 0.8394172690083436, 0.857835825547267, 0.8767043634220298, 0.8959497359587114, 0.9056395781831735, 0.9153294376135723, 0.9249101362920165, 0.9344910292404003, 0.9530118709420836, 0.9704267632055288, 0.986250311068824, 1.0};
*/


real[] reds;
reds.push(.0);
for(int i=1;i<n1-1;++i)
  reds.push(unitrand());
reds.push(1.);

real[] greens;
greens.push(.0);
for(int i=1;i<n2-1;++i)
  greens.push(unitrand());
greens.push(1.);



//regular test with inteleaved reds and greens with gradually squeezing distance
/*
reds.push(.0);
greens.push(.0);
real x=0.;
for(int i=1;i<n1;++i) {
  x += 0.2/i;
  
  if(i%2==0)
    reds.push(x);
  else
    greens.push(x);
}
reds.push(1.);
greens.push(1.);
*/

n1=reds.length;
n2=greens.length;

reds=sort(reds);
greens=sort(greens);

//for(int i=0;i<min(n1,n2);++i)
//  draw((reds[i],0)--(greens[i],1), gray);


//get pairs:



int[][] make_links(int rbegin, int gbegin, int rend, int gend, int dir=1) {
  int[][] linked;
  //write("make_links: ("+string(rbegin)+", "+string(gbegin)+")--("+string(rend)+", "+string(gend)+") dir="+string(dir));
  if(dir>0 && (rbegin>=rend || gbegin>=gend))
    return linked;
  if(dir<0 && (rbegin<=rend || gbegin<=gend))
    return linked;
  
  int candidate[] = {-1,-1};
  bool leading_reds=false;
  int i= rbegin;
  int j= gbegin;
  real distance = 1.;
  bool proceed = true; // helps to run the branch that adds link last time when i and j are already invalid but we still have valid candidate to add. with this boolean we can run the branch in a regular fashion.
  while(true) {
    //write((i,j));
    if(proceed && abs(reds[i] - greens[j]) < distance) {
      distance = abs(reds[i]-greens[j]);
      if(candidate[0]!=-1 && candidate[1]!=-1 && abs(reds[i]-greens[j])<abs(reds[candidate[0]]-greens[candidate[1]])) {
        candidate[0] = i;
        candidate[1] = j;
      }
      if (candidate[0]==-1 && candidate[1]==-1) {
        candidate[0] = i;
        candidate[1] = j;
      }
    } else {
      bool reset=false;
      if(candidate[0]!=-1 && candidate[1]!=-1) {
        //write("candidate: "+string(candidate[0])+", "+string(candidate[1])+" dir="+string(dir));
        int step = dir > 0 ? 1 : -1;
        
        int[][] newlinked = make_links(candidate[0]-step, candidate[1]-step, linked.length == 0 ? rbegin : linked[linked.length-1][0], linked.length == 0 ? gbegin : linked[linked.length-1][1], dir=-dir*2);
        if(newlinked.length != 0) {
          for(int id = 0; id < newlinked.length; ++id) {
            int id_ = dir > 0 ? newlinked.length - 1 - id : id;
            linked.push(newlinked[id_]);
          }
        }
        linked.push(new int[] {candidate[0], candidate[1]});
      }
      candidate[0]=-1;
      candidate[1]=-1;
      distance = 1.;
    }
    if(!proceed)
      break;
    //write("("+string(i)+", "+string(j)+") leading " + (leading_reds?"reds":"greens")+" distance="+string(distance));

    leading_reds = dir * (leading_reds ? 1 : -1) * (reds[i] - greens[j]) > 0 ? !leading_reds : leading_reds;

    int step = dir > 0 ? 1 : -1;
    i=leading_reds ? i+step : i;
    j=!leading_reds ? j+step : j;

    if(i==rend || j==gend)
      proceed=false;

  }
  //write(linked);
  //write("leaving make links");
  return linked;
}

int[][] linked = make_links(0, 0, reds.length, greens.length);

//interval: current + target + bounds
struct interval {
  real[] reds; // to be removed
  real[] greens; // to be added
  real[] added;
  int[] left={0,0};
  int[] right={0,0};
};
interval[] intervals;

//fill intervals
for(int i=0; i<linked.length-1; ++i) {
  interval inter;
  inter.left[0]=linked[i][0];
  inter.left[1]=linked[i][1];
  inter.right[0]=linked[(i+1)%linked.length][0];
  inter.right[1]=linked[(i+1)%linked.length][1];
  for(int j=linked[i][0]+1;j<linked[(i+1)%linked.length][0];++j)
    inter.reds.push((reds[j%reds.length]-reds[inter.left[0]])/(reds[inter.right[0]]-reds[inter.left[0]]));
  for(int j=linked[i][1]+1;j<linked[(i+1)%linked.length][1];++j)
    inter.greens.push((greens[j%greens.length]-greens[inter.left[1]])/(greens[inter.right[1]]-greens[inter.left[1]]));
  intervals.push(inter);
}

interval[][] sections;
sections.push(intervals);

while(true) {
  //make layer
  int removed=0;
  int added=0;
  interval[] newsection;
  for(interval i : sections[sections.length-1]) {
    interval newi;
    newi.left = i.left; 
    newi.right = i.right;
    //loop through current and remove order checkerboard
    real[] newreds;
    for(int j=0;j<i.reds.length;++j) {
      if(j % 2 != 0)
        newreds.push(i.reds[j]);
      else
        removed += 1;
    }
    newi.reds = newreds;
    newi.added = i.added;
    // added+current-target should be monotonic
    if(i.greens.length != 0) {
      real[] ext_interval;
      ext_interval.push(0);
      for(real a : i.added)
        ext_interval.push(a);
      ext_interval.push(1);
      real[] newadded;
      real[] newgreens;
      int j=0;
      int jj=0;
      for(int k=0;k<ext_interval.length-1;++k) {
        int j0 = j;
        int j_to_add = j;
        real median = .5 * (ext_interval[k]+ext_interval[k+1]);
        //write("   "+string(ext_interval[k])+"--"+string(ext_interval[k+1])+", median="+string(median));
        bool found = false;
        for(;j<i.greens.length && i.greens[j] < ext_interval[k+1];++j) {
          //write(i.greens[j]);
          if(!found && i.greens[j] > median) {
            if(j>j0)
              j_to_add = (i.greens[j] - median) < (median - i.greens[j-1]) ? j : j-1;
            else
              j_to_add=j;
            found = true;
          }
        }
        if(!found && j>j0) {
          j_to_add = j-1;
          found = true;
        }
        if(!found) {
          for(;jj<i.added.length && i.added[jj] < ext_interval[k+1]; ++jj)
            newadded.push(i.added[jj]);
          continue;
        }
        //write(j_to_add);
        
        //int jj = 0;
        for(;jj<i.added.length && i.added[jj] < ext_interval[k+1] && i.added[jj]<i.greens[j_to_add]; ++jj)
          newadded.push(i.added[jj]);
        newadded.push(i.greens[j_to_add]);
        for(;jj<i.added.length && i.added[jj] < ext_interval[k+1]; ++jj)
          newadded.push(i.added[jj]);

        added += 1;
        
        for(int jjj=j0;jjj<j;++jjj) {
          if(jjj!=j_to_add)
            newgreens.push(i.greens[jjj]);
        }
      }
      //write("------------");
      newi.added=newadded;
      newi.greens = newgreens;
    }
    newsection.push(newi);
  }
  //write("=======================");
  //write("removed="+string(removed)+", added="+string(added));
  if(removed==0 && added==0)
    break;
  sections.push(newsection);
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
  //write("i0_stride=", i0_stride);
  //write("i1_stride=", i1_stride);
  int i0=0;
  int i1=0;
  int ii1=0;
  int ii0=0;
  for(int i=0;i<intervals0.length;++i) {
    ii1=0;
    ii0=0;
    if(intervals0[i].reds.length != 0) {
      for(; ii1<intervals1[i].reds.length; ++ii1, ii0+=2) {
        int[] q;
        q.push(i0+ii0); q.push(i0+ii0+1); q.push(i0_stride+i1+ii1+1); q.push(i0_stride+i1+ii1);
        indexes.push(q);
        int[] t;
        t.push(i0+ii0+1); t.push(i0+ii0+2); t.push(i0_stride+i1+ii1+1);
        indexes.push(t);
      }
      if(ii0 == intervals0[i].reds.length-1) {
        int[] t;
        t.push(i0+ii0+1); t.push((i0+ii0+2) % i0_stride); t.push(i0_stride+(i1+ii1+1) % i1_stride);
        indexes.push(t);
      }
      int[] q;
      q.push(i0+ii0); q.push((i0+ii0+1) % i0_stride); q.push(i0_stride+(i1+ii1+1) % i1_stride); q.push(i0_stride+(i1+ii1) % i1_stride);
      indexes.push(q);
    }
    if(intervals1[i].added.length != 0) {
      for(; ii0<intervals0[i].added.length; ++ii0, ii1+=2) {
        int[] t;
        t.push(i0+ii0); t.push(i0_stride+i1+ii1+1); t.push(i0_stride+i1+ii1);
        indexes.push(t);
        int[] q;
        q.push(i0+ii0); q.push(i0+ii0+1); q.push(i0_stride+i1+ii1+2); q.push(i0_stride+i1+ii1+1);
        indexes.push(q);
      }
      
      if(ii1 == intervals1[i].added.length-1) {
        int[] t;
        t.push(i0+ii0); t.push(i0_stride+i1+ii1+1); t.push(i0_stride+i1+ii1);
        indexes.push(t);
        int[] q;
        q.push(i0+ii0); q.push((i0+ii0+1) % i0_stride); q.push(i0_stride+(i1+ii1+2) % i1_stride); q.push(i0_stride+i1+ii1+1);
        indexes.push(q);
      } else {
        int[] q;
        q.push(i0+ii0); q.push((i0+ii0+1) % i0_stride); q.push(i0_stride+(i1+ii1+1) % i1_stride); q.push(i0_stride+i1+ii1);
        indexes.push(q);
      }
      
    }
    if(intervals1[i].added.length==0 && intervals0[i].reds.length==0) {
      int[] q;
      q.push(i0+ii0); q.push((i0+ii0+1) % i0_stride); q.push(i0_stride+(i1+ii1+1) % i1_stride); q.push(i0_stride+i1+ii1);
      indexes.push(q);
    }
    i0 += intervals0[i].reds.length + intervals0[i].added.length + 1;
    i1 += intervals1[i].reds.length + intervals1[i].added.length + 1;
  }
  return indexes;
}


void draw_interval_nodes(interval[] intervals, real progress/*0...1*/, real hmax, real shft) {
  for(int id=0; id<intervals.length; ++id) {
    interval i = intervals[id];
    real x0 = reds[i.left[0]]*(1-progress)+greens[i.left[1]]*progress;
    real x1 = reds[i.right[0]]*(1-progress)+greens[i.right[1]]*progress;

    dot((x0,progress*hmax+shft));
    dot((x1,progress*hmax+shft));
    draw((x0+0.01,progress*hmax+shft)--(x1-0.01,progress*hmax+shft));
    for(real v : i.reds)
      dot((v * (x1-x0) + x0,progress*hmax+shft), red);
    for(real v : i.added)
      dot((v * (x1-x0) + x0,progress*hmax+shft), 0.7*blue+0.3*green);
    for(real v : i.greens)
      dot((v * (x1-x0) + x0,progress*hmax+shft-0.003), 0.5*green);
  }
}

real H=0.3;

for(int i=0;i<n1;++i) {
  dot((reds[i],0), red);
  label(rotate(90)*string(i),(reds[i],0),S);
}

for(int i=0;i<n2;++i) {
  dot((greens[i],H), 0.5*green);
  label(rotate(90)*string(i),(greens[i],H),N);
}

for(int[] p : linked)
  draw((reds[p[0]],0)--(greens[p[1]], H), gray);

int layer=0;
for(interval[] intervals : sections) {
  draw_interval_nodes(intervals, layer/(sections.length-1), H-0.01, 0.005);
  layer+=1;
}
