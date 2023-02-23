import bpy
import copy
import bmesh
import mathutils
import copy
import numpy
import math

class interval:
  def __init__(self):
      self.reds=[] #to be removed
      self.greens=[] #to be added
      self.added=[] #actually added
      self.left=[0, 0] #index of left edge, 0-bottom contour, 1-top contour
      self.right=[0, 0] #index of rigth edge, 0-bttom contour, 1-top contour

def connect_intervals_simple(intervals0):
  indexes=[]
  i0_stride=0
  for i in intervals0: 
    i0_stride += 1 + len(i.reds) + len(i.added)
  i0=0
  for i in intervals0:
    for ii in range(0, len(i.reds)):
      q = []
      q.append(i0 + ii)
      q.append((i0 + ii + 1) % i0_stride)
      q.append(i0_stride + (i0 + ii + 1) % i0_stride)
      q.append(i0_stride + i0 + ii)
      indexes.append(q)
    i0+=len(i.reds)

    for ii in range(0, len(i.added)):
      q = []
      q.append(i0 + ii)
      q.append((i0 + ii + 1) % i0_stride)
      q.append(i0_stride + (i0 + ii + 1) % i0_stride)
      q.append(i0_stride + i0 + ii)
      indexes.append(q)
    i0+=len(i.added)

    q = []
    q.append(i0)
    q.append((i0 + 1) % i0_stride)
    q.append(i0_stride + (i0 + 1) % i0_stride)
    q.append(i0_stride + i0)
    indexes.append(q)
    i0+=1
    
  return indexes

def get_params_for_interval(i) :
  res=[.0]
  for  v in i.reds :
    res.append(v)
  for v in i.added :
    res.append(v)
  res.append(1.)
  return res

#returns indexes of quads and tris for pair of intervals of the same length
def connect_intervals(intervals0, intervals1) :
  indexes=[]
  if len(intervals0) != len(intervals1) :
    return indexes
  i0_stride=0
  for i in intervals0 :
    i0_stride += 1 + len(i.reds) + len(i.added)
  i1_stride=0
  for i in intervals1 : 
    i1_stride += 1 + len(i.reds) + len(i.added)
  i0=0
  i1=0
  for i in range(0, len(intervals0)) :
    ii1=0
    ii0=0
    inter0ext = get_params_for_interval(intervals0[i])
    inter1ext = get_params_for_interval(intervals1[i])
    while ii0 < len(inter0ext) - 1 or ii1 < len(inter1ext) - 1 :
      dx11 = 1. if (ii0+1==len(inter0ext) or ii1+1==len(inter1ext)) else abs(inter0ext[ii0+1] - inter1ext[ii1+1])
      dx10 = 1. if ii0+1==len(inter0ext) else abs(inter0ext[ii0+1] - inter1ext[ii1])
      dx01 = 1. if ii1+1==len(inter1ext) else abs(inter0ext[ii0] - inter1ext[ii1+1])
      if dx11 <= dx10 and dx11 <= dx01 :
        q = []
        q.append(i0+ii0)
        q.append((i0+ii0+1) % i0_stride)
        q.append(i0_stride+(i1+ii1+1) % i1_stride)
        q.append(i0_stride+i1+ii1)
        indexes.append(q)
        ii0 += 1
        ii1 += 1
      if dx10 <= dx11 and dx10 <= dx01 :
        t = []
        t.append(i0+ii0)
        t.append((i0+ii0+1) % i0_stride)
        t.append(i0_stride+(i1+ii1) % i1_stride)
        indexes.append(t)
        ii0 += 1
      if dx01 <= dx11 and dx01 <= dx10 :
        t = []
        t.append((i0+ii0) % i0_stride)
        t.append(i0_stride+(i1+ii1+1) % i1_stride)
        t.append(i0_stride+i1+ii1)
        indexes.append(t)
        ii1 += 1
    i0 += len(inter0ext) - 1
    i1 += len(inter1ext) - 1
  return indexes

def make_sections(reds, greens):
  #a and b are extended contours parametrization i.e. all points are included [0...1], point with parameter 1 will be treated as 0 point duplicate

    def make_links(rbegin, gbegin, rend, gend, dir=1) :
      linked=[]
      #print("make_links: ("+str(rbegin)+", "+str(gbegin)+")--("+str(rend)+", "+str(gend)+") dir="+string(dir));
      if dir>0 and (rbegin>=rend or gbegin>=gend) :
        return linked
      if dir<0 and (rbegin<=rend or gbegin<=gend) :
        return linked
      
      candidate=[-1,-1]
      leading_reds=False
      i = rbegin
      j = gbegin
      distance = 1.
      proceed = True # helps to run the branch that adds link last time when i and j are already invalid but we still have valid candidate to add. with this boolean we can run the branch in a regular fashion.
      while True :
        #print("("+str(i)+", "+str(j)+")")
        if proceed and abs(reds[i] - greens[j]) < distance :
          distance = abs(reds[i]-greens[j])
          if candidate[0] != -1 and candidate[1] != -1 and abs(reds[i]-greens[j])<abs(reds[candidate[0]]-greens[candidate[1]]) :
            candidate[0] = i
            candidate[1] = j
          if candidate[0] == -1 and candidate[1] == -1 :
            candidate[0] = i
            candidate[1] = j
        else :
          reset=False
          if candidate[0] != -1 and candidate[1] != -1 :
            #print("candidate: "+str(candidate[0])+", "+str(candidate[1])+" dir="+str(dir))
            step = 1 if dir > 0 else -1
            newlinked = make_links(candidate[0] - step, candidate[1]-step, rbegin if len(linked) == 0 else linked[len(linked)-1][0], gbegin if len(linked) == 0 else linked[len(linked) - 1][1], dir = -dir * 2)
            if len(newlinked) != 0 :
              for id in range(0, len(newlinked)) :
                id_ = (len(newlinked) - 1 - id) if dir > 0 else id
                linked.append(newlinked[id_])
            linked.append([candidate[0], candidate[1]]);
          candidate[0] = -1
          candidate[1] = -1
          distance = 1.
        if not proceed :
          break
        #print("("+str(i)+", "+str(j)+") leading " + ("reds" if leading_reds==True else "greens")+" distance="+str(distance))

        leading_reds = (not leading_reds) if (dir * (1 if leading_reds==True else -1) * (reds[i] - greens[j]) > 0) else leading_reds

        step = 1 if dir > 0 else -1
        i = i+step if leading_reds==True else i;
        j = j+step if leading_reds==False else j

        if i==rend or j==gend:
          proceed=False

      #print(linked)
      #print("leaving make links")
      return linked
  
    linked = make_links(0, 0, len(reds), len(greens))

    #print(linked)

    intervals = []

    for i in range(0, len(linked)-1) :
      inter = interval()
      inter.left[0] = linked[i][0]
      inter.left[1] = linked[i][1]
      inter.right[0] = linked[(i+1) % len(linked)][0]
      inter.right[1] = linked[(i+1) % len(linked)][1]
      for j in range(linked[i][0]+1, linked[(i+1) % len(linked)][0]) :
        inter.reds.append((reds[j % len(reds)]-reds[inter.left[0]])/(reds[inter.right[0]]-reds[inter.left[0]]))
      for j in range(linked[i][1]+1, linked[(i+1) % len(linked)][1]) :
        inter.greens.append((greens[j % len(greens)]-greens[inter.left[1]])/(greens[inter.right[1]]-greens[inter.left[1]]))
      intervals.append(inter)


    sections = []
    sections.append(intervals)
    
    while True :
      #make layer
      removed=0
      added=0
      newsection = []
      for i in sections[len(sections) - 1] :
        #print("***")
        #print(i.added)
        newi = interval()
        newi.left = i.left 
        newi.right = i.right
        #loop through current and remove order checkerboard
        newreds=[]
        for j in range(0, len(i.reds)) :
          if j % 2 != 0 :
            newreds.append(i.reds[j])
          else :
            removed += 1
        newi.reds = newreds
        newi.added = i.added
        #added+current-target should be monotonic
        if len(i.greens) != 0:
          ext_interval = []
          ext_interval.append(.0)
          for a in i.added :
            ext_interval.append(a)
          ext_interval.append(1.)
          newadded = []
          newgreens = []
          j = 0
          jj = 0
          for k in range(0, len(ext_interval) - 1) :
            j0 = j
            j_to_add = j
            median = .5 * (ext_interval[k]+ext_interval[k+1])
            #print("   "+str(ext_interval[k])+"--"+str(ext_interval[k+1])+", median="+str(median))
            found = False;
            while j < len(i.greens) and i.greens[j] < ext_interval[k+1] :
              #print("i.greens["+str(j)+"]="+str(i.greens[j]))
              if (not found) and i.greens[j] > median :
                if j > j0 :
                  j_to_add = j if (i.greens[j] - median) < (median - i.greens[j-1]) else j-1;
                else :
                  j_to_add = j
                found = True;
              j += 1
            if not found and j > j0 :
              j_to_add = j-1
              found = True;
            if not found:
              #print("not found")
              while jj < len(i.added) and i.added[jj] < ext_interval[k+1] :
                newadded.append(i.added[jj])
                jj += 1
              continue
            #else:
            #  print("found")
            #print("j_to_add="+str(j_to_add)+" value=" + str(i.greens[j_to_add]))
            
            #jj = 0
            while jj < len(i.added) and i.added[jj] < ext_interval[k+1] and i.added[jj] < i.greens[j_to_add] : 
              newadded.append(i.added[jj])
              jj += 1
            newadded.append(i.greens[j_to_add])
            while jj < len(i.added) and i.added[jj] < ext_interval[k+1] :
              newadded.append(i.added[jj])
              jj += 1

            added += 1
            
            for jjj in range(j0, j) :
              if jjj != j_to_add :
                newgreens.append(i.greens[jjj])
          #print("------------")
          newi.added = newadded
          #print(newi.added)
          newi.greens = newgreens
        newsection.append(newi)
      #print("=======================")
      #print("removed="+str(removed)+", added="+str(added))
      sections.append(newsection)
      if removed == 0 and added == 0 :
        break
    return sections
  
#SEAM METHODS
def sortListFirstItemFunc(a):
  return a[0]

#subdivides s and returns new subdevided s
#s is a parametrization
def subdivide(s, n_desired, accurate=True):
  s_sub = s

  while len(s_sub) < n_desired :
    n = len(s_sub)
    if n >= n_desired:
      return s_sub

    seglenghts = []

    for pid in range(0,n-1):
      seglenghts.append([s_sub[pid+1] - s_sub[pid], pid])
    
    seglenghts.sort(key=sortListFirstItemFunc, reverse=True)

    ids=[]
    for i in range(0,len(seglenghts)) :
      ids.append(0)
    i=0
    while n + i < n_desired and i < len(seglenghts) :
      ids[seglenghts[i][1]] = 1
      i+=1
      if accurate==True : #subdivide only biggest interval on each iteration
        break
    new_s = []
    for pid in range(0, n - 1) :
      new_s.append(copy.deepcopy(s_sub[pid]))
      if ids[pid] == 1 :
        new_s.append((s_sub[pid+1] + s_sub[pid]) * .5)
    new_s.append(copy.deepcopy(s_sub[n-1]))
    s_sub = new_s
  
  return s_sub

def get_ordered_verts(obj):
  cont = []
  if obj.data.bl_rna == bpy.types.Curve.bl_rna :
    #print("obj is spline")
    spl=obj.data.splines[0]
    minus_one = -1 if spl.use_cyclic_u == False else 0
    for sid in range(0, len(spl.bezier_points) + minus_one):
      sid_next = (sid + 1 ) % len(spl.bezier_points)
      seg_points = mathutils.geometry.interpolate_bezier(obj.matrix_world * spl.bezier_points[sid].co, obj.matrix_world * spl.bezier_points[sid].handle_right, obj.matrix_world * spl.bezier_points[sid_next].handle_left, obj.matrix_world * spl.bezier_points[sid_next].co, spl.resolution_u)     
      for pid in range(0,len(seg_points)+(-1 if sid != len(spl.bezier_points)-2 else 0)):
        cont.append(seg_points[pid])
    
  else :
    #print("obj is mesh")
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    for v0 in bm.verts:
      vp=v0
      v=v0.link_edges[0].other_vert(v0)
      cont.append(obj.matrix_world * vp.co)
      while v != v0:
        cont.append(obj.matrix_world * v.co)
        nextV = v.link_edges[0].other_vert(v)
        if len(v.link_edges) == 1:
          break
        if nextV == vp:
          nextV = v.link_edges[1].other_vert(v)
        vp = v
        v = nextV
      break
    bm.free()
  return cont

def sortListFirstItemFunc(a):
  return a[0]

def loopindex(sz, i):
  return i % sz

def loopvert(obj, i):
  return obj[loopindex(len(obj),i)]

def closest_vert(s, c):
  res = [-1,-1]
  minr2 = 10000
  for i in range(0,len(s)):
    coord = s[i]
    r2 = 100000
    pre_res_1 = -1
    for j in range(0,len(c)):
      v = c[j]
      r = coord - v
      curr2 = r.length_squared
      if curr2 < r2:
        pre_res_1 = j
        r2 = curr2
    if r2 < minr2:
      res[0] = i
      res[1] = pre_res_1
      minr2 = r2
  return res

def get_orths(seam, debug=False):
  #edge orths:
  ex0s = []
  ey0s = []
  ez0s = []
  ez0s.append((seam[1]-seam[0]).normalized())
  ey0s.append(ez0s[0].orthogonal().normalized())
  ex0s.append(ez0s[0].cross(ey0s[0]))  
  for ei in range(1,len(seam)-1):
    ez0 = ez0s[-1]
    ez0new = (seam[ei+1]-seam[ei]).normalized()
    n = ez0.cross(ez0new)
    if n.length_squared<0.0001*ez0new.length_squared : #we cannot determine orth for straight path so we use previous values
      ez0s.append(ez0s[-1])
      ey0s.append(ey0s[-1])
      ex0s.append(ex0s[-1])
      continue

    n.normalize()
    nperp = n.cross(ez0)
    #nnew = n //by method
    nperpnew = n.cross(ez0new)
    ey0 = ey0s[-1]
    ey0znperp = ey0 - n.dot(ey0)*n
    ey0znperpew = nperp.dot(ey0znperp)* nperpnew
    ey0new = ey0znperpew + n.dot(ey0)*n
    ex0new = ez0new.cross(ey0new)
    ex0s.append(ex0new)
    ey0s.append(ey0new)
    ez0s.append(ez0new)
  
  if debug:
      mesh = bpy.data.meshes.new("mesh")
      target = bpy.data.objects.new("seam local coords e", mesh)
      bpy.context.scene.objects.link(target)
      bm = bmesh.new()
      for i in range(0,len(ez0s)):
        origin=0.5*(seam[i]+seam[i+1])
        factor=0.1
        v0 = bm.verts.new(origin)
        v1 = bm.verts.new(origin+factor*ex0s[i])
        v2 = bm.verts.new(origin+factor*ey0s[i])
        v3 = bm.verts.new(origin+factor*ez0s[i])
        bm.edges.new((v0,v1))
        bm.edges.new((v0,v2))
        bm.edges.new((v0,v3))
      bm.to_mesh(mesh)  
      bm.free()
  #interpolate orths from edges(edge centers) to seam vertices
  x0s = []
  y0s = []
  z0s = []
  x0s.append(ex0s[0])
  y0s.append(ey0s[0])
  z0s.append(ez0s[0])
  n=len(seam)
  for i in range(0,n-2):
    z0=0.5*(ez0s[i]+ez0s[i+1])
    z0.normalize()
    x0=0.5*(ex0s[i]+ex0s[i+1])
    y0 = x0.cross(z0)
    y0.normalize()
    x0 = z0.cross(y0)
    x0s.append(x0)
    y0s.append(y0)
    z0s.append(z0)
  x0s.append(ex0s[n-2])
  y0s.append(ey0s[n-2])
  z0s.append(ez0s[n-2])

  if debug:
      mesh = bpy.data.meshes.new("mesh")
      target = bpy.data.objects.new("seam local coords v", mesh)
      bpy.context.scene.objects.link(target)
      bm = bmesh.new()
      for i in range(0,len(x0s)):
        origin=seam[i]
        factor=0.1
        v0 = bm.verts.new(origin)
        v1 = bm.verts.new(origin+factor*x0s[i])
        v2 = bm.verts.new(origin+factor*y0s[i])
        v3 = bm.verts.new(origin+factor*z0s[i])
        bm.edges.new((v0,v1))
        bm.edges.new((v0,v2))
        bm.edges.new((v0,v3))
      bm.to_mesh(mesh)  
      bm.free()    
  return [x0s,y0s,z0s]

def get_S(c, i=0):
  s=mathutils.Vector((0,0,0))
  for j in range(1,len(c)-1):
    a=loopvert(c,i+j)-loopvert(c,i)
    b=loopvert(c,i+j+1)-loopvert(c,i)
    s+=b.cross(a)
  return s

def get_anchors(contours, seam, debug=False):
    anchors = []
    for ci in range(0,len(contours)):
      sici0 = closest_vert(seam,contours[ci])
      sici0.append(ci)
      anchors.append(sici0)

    if debug:
      print(anchors)
    
    anchors.sort(key=sortListFirstItemFunc)
    return anchors

def  fix_contours_orientation(anchors, contours, seam): #CAUTION: changes contours
    #fix contours orientation
    for a in anchors:
      dir=seam[a[0]+1 if a[0]+1 < len(seam) else a[0]]-seam[a[0]-1 if a[0]-1 >= 0 else a[0]]
      cont = contours[a[2]]
      if get_S(cont,a[1]).dot(dir) < 0:
        newcont = copy.deepcopy(cont)
        for i in range(0,len(cont)):
          newcont[loopindex(len(newcont),a[1]-i)] = loopvert(cont, a[1]+i)
        contours[a[2]] = newcont

def spline(xs, ys, x, derivative0=None, derivative1=None, linear=False):
  res = [.0,.0,.0]

  if len(xs)<4 or len(ys)<4 or x<xs[1] or x>xs[2]:
    return ys[0]

  if linear:
    matrix = numpy.array([[xs[1],1], [xs[2],1]])
    b = numpy.array([ys[1], ys[2]])
    koeffs = numpy.linalg.solve(matrix,b).tolist()
    for i in range(0,3):
      res[i] = x*koeffs[0][i]+koeffs[1][i]
    return mathutils.Vector((res[0],res[1],res[2]))

  have_deriv0 = False
  deriv0 = ys[0]
  if derivative0 == None :
   if xs[0]!=xs[1]:
     deriv0 = (1/(xs[2]-xs[0]))*(ys[2]-ys[0])
     have_deriv0 = True
  else :
    have_deriv0 = True
    deriv0 = derivative0
  
  deriv1 = ys[0]
  have_deriv1 = False
  if derivative1 == None :
    if xs[2]!=xs[3]:
      deriv1 = (1/(xs[3]-xs[1]))*(ys[3]-ys[1])
      have_deriv1 = True
  else :
    have_deriv1 = True
    deriv1 = derivative1

  if have_deriv0 and have_deriv1:
    matrix = numpy.array([[xs[1]**3,xs[1]**2,xs[1],1], [xs[2]**3,xs[2]**2,xs[2],1], [3*xs[1]**2,2*xs[1],1,0], [3*xs[2]**2,2*xs[2],1,0]])
    b = numpy.array([ys[1], ys[2], deriv0, deriv1])
    koeffs = numpy.linalg.solve(matrix,b).tolist()
    for i in range(0,3):
      res[i] = x**3*koeffs[0][i]+x**2*koeffs[1][i]+x*koeffs[2][i]+koeffs[3][i]
  if not have_deriv0 and have_deriv1:
    matrix = numpy.array([[xs[1]**2,xs[1],1], [xs[2]**2,xs[2],1], [2*xs[2],1,0]])
    b=numpy.array([ys[1], ys[2], deriv1])
    koeffs = numpy.linalg.solve(matrix,b).tolist()
    for i in range(0,3):
      res[i] = x**2*koeffs[0][i]+x*koeffs[1][i]+koeffs[2][i]
  if have_deriv0 and not have_deriv1:
    matrix = numpy.array([[xs[1]**2,xs[1],1], [xs[2]**2,xs[2],1], [2*xs[1],1,0]])
    b=numpy.array([ys[1], ys[2], deriv0])
    koeffs = numpy.linalg.solve(matrix,b).tolist()
    for i in range(0,3):
      res[i] = x**2*koeffs[0][i]+x*koeffs[1][i]+koeffs[2][i]
  if not have_deriv0 and not have_deriv1:
    matrix = numpy.array([[xs[1],1], [xs[2],1]])
    b = numpy.array([ys[1], ys[2]])
    koeffs = numpy.linalg.solve(matrix,b).tolist()
    for i in range(0,3):
      res[i] = x*koeffs[0][i]+koeffs[1][i]
  return mathutils.Vector((res[0],res[1],res[2]))

def get_normal(c) :
    n = mathutils.Vector((0.,0.,0.))
    for v in c:
      n += (c[1] - c[0]).cross(v - c[0])
    return n

def ShowMessageBox(message = "", title = "Error", icon = 'ERROR'):
    def draw(self, context):
        self.layout.label(text=message)
    bpy.context.window_manager.popup_menu(draw, title = title, icon = icon)

def make_tube(contours, seam, alpha0_deg=None, in_local_space=True, alpha1_deg=None, linear=False) :
    if alpha0_deg != None:
      alpha0 = alpha0_deg / 180 * math.pi
    else:
      alpha0 = None
    if alpha1_deg != None:
      alpha1 = alpha1_deg / 180 * math.pi
    else:
      alpha1 = None
    nc=len(contours[0])
    allgood=True
    for c in contours:
      if len(c)!=nc :
        allgood=False
    if allgood == False:
      ShowMessageBox('Contours must be of the same size')
      return

    local_space_chain=get_orths(seam)
    
    def to_local_space(v, i):
      if in_local_space == True :
        return mathutils.Vector((v.dot(local_space_chain[0][i]),v.dot(local_space_chain[1][i]),v.dot(local_space_chain[2][i])))
      else:
        return v
      
    def to_global_space(v, i):
      if in_local_space == True :
        return (v.x*local_space_chain[0][i]+v.y*local_space_chain[1][i]+v.z*local_space_chain[2][i])
      else:
        return v
    
    seam_path_len = [0.]
    for vi in range(1,len(seam)):
      seam_path_len.append(seam_path_len[-1] + (seam[vi]-seam[vi-1]).length)

    anchors = get_anchors(contours, seam, debug=False)

    fix_contours_orientation(anchors, contours, seam)

    bm = bmesh.new()

    n0 = get_normal(contours[anchors[0][2]])
    n_1 = get_normal(contours[anchors[-1][2]])

    sect_count=1
    contours_length = len(contours[0]) #must be the same for all the contours
    allverts = []
    for ai in range(0,len(anchors)-1):
      ai_1=ai-1 if ai>0 else 0 
      ai0=ai
      ai1=ai+1
      ai2=ai+2 if ai+2<len(anchors) else ai+1
      
      i_1 = anchors[ai_1][0]
      i0 = anchors[ai0][0]
      i1 = anchors[ai1][0]
      i2 = anchors[ai2][0]
      x_1=seam_path_len[i_1]
      x0=seam_path_len[i0]
      x1=seam_path_len[i1]
      x2=seam_path_len[i2]
      for i in range(i0, i1):
        j0_1=anchors[ai_1][1]
        j00=anchors[ai0][1]
        j01=anchors[ai1][1]
        j02=anchors[ai2][1]
        for j in range(0,contours_length):
              
          c_1vloc = loopvert(contours[anchors[ai_1][2]],j+j0_1) - seam[i_1]
          c0vloc = loopvert(contours[anchors[ai0][2]],j+j00) - seam[i0]
          c1vloc = loopvert(contours[anchors[ai1][2]],j+j01) - seam[i1]
          c2vloc = loopvert(contours[anchors[ai2][2]],j+j02) - seam[i2]
          
          c_1vloc = to_local_space(c_1vloc, i_1)
          c0vloc = to_local_space(c0vloc, i0)
          c1vloc = to_local_space(c1vloc, i1)
          c2vloc = to_local_space(c2vloc, i2)

          derivative0 = None
          if ai==0 and alpha0 != None:
            v1 = to_local_space(loopvert(contours[anchors[ai0][2]],j+j00-1)-seam[i0],i0) - c0vloc
            v2 = to_local_space(loopvert(contours[anchors[ai0][2]],j+j00+1)-seam[i0],i0) - c0vloc
            #n=v1.cross(v2)
            n0_loc = to_local_space(n0, i0)
            tau=v1-v2
            r=n0_loc.cross(tau)
            derivative0=(1./r.length) * math.tan(alpha0) * r
          derivative1 = None
          if ai==len(anchors)-2 and alpha1 != None:
            v1 = to_local_space(loopvert(contours[anchors[ai1][2]],j+j01-1)-seam[i1],i1) - c1vloc
            v2 = to_local_space(loopvert(contours[anchors[ai1][2]],j+j01+1)-seam[i1],i1) - c1vloc
            #n=v1.cross(v2)
            n_1_loc = to_local_space(n_1, i1)
            tau=v1-v2
            r=n_1_loc.cross(tau)
            derivative1=(1./r.length) * math.tan(alpha1) * r
                    
          x=seam_path_len[i]
          #pt = (1/(x1-x0))*((x1-x)*c0vloc+(x-x0)*c1vloc)
          pt = spline([x_1,x0,x1,x2], [c_1vloc,c0vloc,c1vloc,c2vloc], x, derivative0 = derivative0, derivative1 = derivative1, linear=linear)

          pt = to_global_space(pt, i)
          allverts.append(bm.verts.new(pt+seam[i]))
        sect_count+=1;
    for j in range(0,contours_length):
      pt = loopvert(contours[anchors[-1][2]],j+anchors[-1][1])
      allverts.append(bm.verts.new(pt))

    for i in range(0, sect_count-1):
      for j in range(0, contours_length):
        v0 = allverts[j+i*contours_length]
        v1 = allverts[((j+1)%contours_length)+i*contours_length]
        v2 = allverts[((j+1)%contours_length)+(i+1)*contours_length]
        v3 = allverts[j+(i+1)*contours_length]
        bm.faces.new((v0,v1,v2,v3))

    mesh = bpy.data.meshes.new("mesh")
    target = bpy.data.objects.new("multi-contour_sweep", mesh)
    bpy.context.scene.objects.link(target)
      
    bm.to_mesh(mesh)  
    bm.free()

def get_contour_parametrization(contour, id0):
  contour_parametrization=[.0]
  total=.0
  for i in range(0, len(contour)):
    v0 = loopvert(contour, id0+i)
    v1 = loopvert(contour, id0+i+1)
    total += (v1 - v0).length

  for i in range(0, len(contour)):
    v0 = loopvert(contour, id0+i)
    v1 = loopvert(contour, id0+i+1)
    contour_parametrization.append(contour_parametrization[-1] + (v1 - v0).length)
  for i in range(0, len(contour_parametrization)):       
    contour_parametrization[i] /= total
  return contour_parametrization

def get_contour_vertex_by_parameter(contour, id0, parameter):
  parametrization = get_contour_parametrization(contour, id0)
  for i in range(0, len(parametrization)):
    if parametrization[i] >= parameter:
      if i == 0:
        return loopvert(contour, 0+id0)
      else:
        p = (parameter - parametrization[i-1]) / (parametrization[i] - parametrization[i-1])
        return (1. - p) * loopvert(contour, i-1+id0) + p * loopvert(contour, i+id0)
  return mathutils.Vector(0., 0.)


def make_t_or_q_face(ind, allverts, vert_id_base, bm) :
  good = True;
  if len(ind) == 3 :
    if vert_id_base + ind[0] >= len(allverts) or vert_id_base + ind[1] >= len(allverts) or vert_id_base + ind[2] >= len(allverts) :
      good=False
    else:
     v0 = allverts[vert_id_base + ind[0]]
     v1 = allverts[vert_id_base + ind[1]]
     v2 = allverts[vert_id_base + ind[2]]
     bm.faces.new((v0, v1, v2))
  if len(ind) == 4 :
    if vert_id_base + ind[0] >= len(allverts) or vert_id_base + ind[1] >= len(allverts) or vert_id_base + ind[2] >= len(allverts) or vert_id_base + ind[3] >= len(allverts) :
      good=False
    else:
     v0 = allverts[vert_id_base + ind[0]]
     v1 = allverts[vert_id_base + ind[1]]
     v2 = allverts[vert_id_base + ind[2]]
     v3 = allverts[vert_id_base + ind[3]]
     bm.faces.new((v0, v1, v2, v3))
  return good

def make_tube_experimental(contours, seam, in_local_space=True, alpha0_deg=None, alpha1_deg=None, linear=False) :
    if alpha0_deg != None:
      alpha0 = alpha0_deg / 180 * math.pi
    else:
      alpha0 = None
    if alpha1_deg != None:
      alpha1 = alpha1_deg / 180 * math.pi
    else:
      alpha1 = None

    local_space_chain=get_orths(seam)
    
    def to_local_space(v, coord_system):
      if in_local_space == True :
        return mathutils.Vector((v.dot(coord_system[0]),v.dot(coord_system[1]),v.dot(coord_system[2])))
      else :
        return v
    
    def to_global_space(v, coord_system):
      if in_local_space == True :
        return (v.x*coord_system[0]+v.y*coord_system[1]+v.z*coord_system[2])
      else :
        return v
    
    seam_path_len = [0.]
    for vi in range(1,len(seam)):
      seam_path_len.append(seam_path_len[-1] + (seam[vi]-seam[vi-1]).length)

    anchors = get_anchors(contours, seam, debug=False)

    fix_contours_orientation(anchors, contours, seam)

    bm = bmesh.new()

    allverts = []
    vert_id_base = 0

    n0 = get_normal(contours[anchors[0][2]])
    n_1 = get_normal(contours[anchors[-1][2]])
    
    for ai in range(0,len(anchors)-1):
      ai_1=ai-1 if ai>0 else 0 
      ai0=ai
      ai1=ai+1
      ai2=ai+2 if ai+2<len(anchors) else ai+1

      reds = get_contour_parametrization(contours[anchors[ai0][2]], anchors[ai0][1])
      greens = get_contour_parametrization(contours[anchors[ai1][2]], anchors[ai1][1])
      
      #red is starting anchor contour within contours pair
      #green is end anchor contour within contour pairs
      sections = make_sections(reds, greens)
      
      #number of sections is equal to number of transitions of contour vertices
      #each section consist of intervals. every section has the same number of intervals
      #within given pair of anchor contours
      #interval is a complex structure introduced above
      sn=len(sections)
      si_1=anchors[ai_1][0]
      si0=anchors[ai0][0]
      si1=anchors[ai1][0]
      si2 = anchors[ai2][0]
      z_1=seam_path_len[si_1]
      z0=seam_path_len[si0]
      z1=seam_path_len[si1]
      z2=seam_path_len[si2]

      #local seam is part of seam for [si0,si1] interval
      #it is used to subdivide seam if needed i.e. if number of sections is greater than si1-si0
      seam_param = [0]
      for si in range(si0, si1) :
        seam_param.append(seam_param[-1] + (seam[si+1] - seam[si]).length)
    
      fine_seam_param = subdivide(seam_param, sn)
      
      local_seam = []
      fsi = 0
      segment_local_space_chain = []
      segment_local_space_chain.append([])
      segment_local_space_chain.append([])
      segment_local_space_chain.append([])
      for si in range(si0,si1) :
        while seam_param[si+1-si0] != fine_seam_param[fsi] :
          p=(fine_seam_param[fsi]-seam_param[si-si0])/(seam_param[si+1-si0]-seam_param[si-si0])
          #print("si="+str(si)+", fsi="+str(fsi)+", p="+str(p))
          local_seam.append((1.-p) * seam[si] + p * seam[si+1])
          for ci in range(0,3) :
            segment_local_space_chain[ci].append((1.-p) * local_space_chain[ci][si] + p * local_space_chain[ci][si+1])
          fsi+=1
      local_seam.append(seam[si1])
      for ci in range(0,3) :
        segment_local_space_chain[ci].append(local_space_chain[ci][si1])
          
      local_seam_path_len = [seam_path_len[si0]]
      si0=0
      si1=len(local_seam)-1
      for si in range(si0,si1) :
        local_seam_path_len.append(local_seam_path_len[-1]+(local_seam[si+1]-local_seam[si]).length)    
            
      sections_parametrization = []
      #parametrization is a vector parameter value for each vertex of section
      #sections_parametrization is a collection of sections parametrization
      #at this point we attach sections to certain z value (along given seam)
      #seam path within [si0, si1] should have >= vertices than number of sections
      for id in range(0, len(sections)):
        section = sections[id]
        parametrization = []
        z = (local_seam_path_len[round(id / (sn - 1) * (si1 - si0))] - z0) / (z1 - z0)
        #NOTE that if si1-si0 > sn-1 then z(id) is never equal to z(id+1)
        x1 = 0
        for inter in section:
          x0 = reds[inter.left[0]] * (1. - z) + greens[inter.left[1]] * z
          x1 = reds[inter.right[0]] * (1. - z) + greens[inter.right[1]] * z
          parametrization.append(x0)
          for p in inter.added:
            parametrization.append(x0 + p * (x1 - x0))
          for p in inter.reds:
            parametrization.append(x0 + p * (x1 - x0))
        parametrization.append(x1)
        sections_parametrization.append(parametrization)

      prev_id=-1
      indexes=[]
      for si in range(si0, si1+1):
        id = int((si-si0)/(si1+1-si0)*sn)
        if id >= sn :
          continue
        z = local_seam_path_len[si-si0] #(local_seam_path_len[si-si0] - z0) / (z1 - z0)
        if prev_id != -1:
          if prev_id != id:
            indexes = connect_intervals(sections[prev_id], sections[id])
          else:
            indexes = connect_intervals_simple(sections[id])
        
        if ai==0 or prev_id != -1 :
          vert_id_base_ = len(allverts)
          npar = len(sections_parametrization[id])
          for pid in range(0, npar):
            if pid == len(sections_parametrization[id]) - 1 :
              continue
            p = sections_parametrization[id][pid]
            v_1 = get_contour_vertex_by_parameter(contours[anchors[ai_1][2]], anchors[ai_1][1],p) - seam[si_1]
            v0 = get_contour_vertex_by_parameter(contours[anchors[ai0][2]], anchors[ai0][1],p) - local_seam[0]
            v1 = get_contour_vertex_by_parameter(contours[anchors[ai1][2]], anchors[ai1][1],p) - local_seam[si1-si0]
            v2 = get_contour_vertex_by_parameter(contours[anchors[ai2][2]], anchors[ai2][1],p) - seam[si2]

            v_1_local_space = [local_space_chain[0][si_1], local_space_chain[1][si_1], local_space_chain[2][si_1]]
            v0_local_space = [segment_local_space_chain[0][0], segment_local_space_chain[1][0], segment_local_space_chain[2][0]]
            v1_local_space = [segment_local_space_chain[0][si1-si0], segment_local_space_chain[1][si1-si0], segment_local_space_chain[2][si1-si0]] 
            v2_local_space = [local_space_chain[0][si2], local_space_chain[1][si2], local_space_chain[2][si2]] 
            
            v_1_loc = to_local_space(v_1, v_1_local_space)
            v0_loc = to_local_space(v0, v0_local_space)
            v1_loc = to_local_space(v1, v1_local_space)
            v2_loc = to_local_space(v2, v2_local_space)

            derivative0 = None
            if ai==0 and alpha0 != None:
              v1 = to_local_space(get_contour_vertex_by_parameter(contours[anchors[ai0][2]], anchors[ai0][1], sections_parametrization[id][(pid-1+npar) % npar]) - local_seam[0], v0_local_space) - v0_loc
              v2 = to_local_space(get_contour_vertex_by_parameter(contours[anchors[ai0][2]], anchors[ai0][1], sections_parametrization[id][(pid+1) % npar]) - local_seam[0], v0_local_space) - v0_loc
              n0_loc = to_local_space(n0, v0_local_space)
              tau=v1-v2
              r=n0_loc.cross(tau)
              derivative0=(1./r.length) * math.tan(alpha0) * r
            derivative1 = None
            if ai==len(anchors)-2 and alpha1 != None:
              v1 = to_local_space(get_contour_vertex_by_parameter(contours[anchors[ai1][2]], anchors[ai1][1], sections_parametrization[id][(pid-1+npar) % npar]) - local_seam[si1-si0], v1_local_space) - v1_loc
              v2 = to_local_space(get_contour_vertex_by_parameter(contours[anchors[ai1][2]], anchors[ai1][1], sections_parametrization[id][(pid+1) % npar]) - local_seam[si1-si0],v1_local_space) - v1_loc
              n_1_loc = to_local_space(n_1, v1_local_space)
              tau=v1-v2
              r=n_1_loc.cross(tau)
              derivative1=(1./r.length) * math.tan(alpha1) * r

            
            v_loc = spline([z_1, z0, z1, z2], [v_1_loc, v0_loc, v1_loc, v2_loc], z, derivative0 = derivative0, derivative1 = derivative1, linear=linear)
            
            v = to_global_space(v_loc, [segment_local_space_chain[0][si], segment_local_space_chain[1][si], segment_local_space_chain[2][si]])
            allverts.append(bm.verts.new(v + local_seam[si-si0]))
        if prev_id != -1:
            good = True
            for ind in indexes:
              good &= make_t_or_q_face(ind, allverts, vert_id_base, bm)
            if good==False :
              print("STOPPPED!")
              break
        vert_id_base = vert_id_base_
        prev_id = id
      #break    
    tube_mesh = bpy.data.meshes.new("mesh")
    target = bpy.data.objects.new("multi-contour_sweep_exp", tube_mesh)
    bpy.context.scene.objects.link(target)    
    bm.to_mesh(tube_mesh)  
    bm.free()



bl_info = {
    "name": "Multi-contour extrusion",
    "category": "Object",
}

class ObjectCursorArray(bpy.types.Operator):
    """Multi-contour extrusion"""
    bl_idname = "object.multi_cont_extrusion"
    bl_label = "Multi-contour extrusion"
    bl_options = {'REGISTER', 'UNDO'}

    exp = bpy.props.BoolProperty(name="experimental", description="extrusion for contours with different number of vertices", default=False)
    lin = bpy.props.BoolProperty(name="linear", description="linear/spline connection between contours", default=False)
    loc = bpy.props.BoolProperty(name="local", description="connect contours in local space", default=True)
    alpha0 = bpy.props.FloatProperty(name="alpha0", description="initial angle", default = 0)
    alpha1 = bpy.props.FloatProperty(name="alpha1", description="final angle", default = 0)

    def execute(self, context):
      conts = []
 
      for i in bpy.context.selected_objects:
        if i != bpy.context.active_object:
          conts.append(get_ordered_verts(i))

      s = get_ordered_verts(bpy.context.active_object)
      if self.exp == True :
        make_tube_experimental(conts, s, in_local_space=self.loc, alpha0_deg=self.alpha0, alpha1_deg=self.alpha1, linear=self.lin)
      else:
        make_tube(conts, s, in_local_space=self.loc, alpha0_deg=self.alpha0, alpha1_deg=self.alpha1, linear=self.lin)   

      return {'FINISHED'}


def menu_func(self, context):
    self.layout.operator(ObjectCursorArray.bl_idname)

# store keymaps here to access after registration
addon_keymaps = []


def register():
    bpy.utils.register_class(ObjectCursorArray)
    bpy.types.VIEW3D_MT_object.append(menu_func)

    # handle the keymap
    wm = bpy.context.window_manager
    # Note that in background mode (no GUI available), keyconfigs are not available either,
    # so we have to check this to avoid nasty errors in background case.
    kc = wm.keyconfigs.addon
    if kc:
        km = wm.keyconfigs.addon.keymaps.new(name='Object Mode', space_type='EMPTY')
        kmi = km.keymap_items.new(ObjectCursorArray.bl_idname, 'SPACE', 'PRESS', ctrl=True, shift=True)
        #kmi.properties.total = 4
        addon_keymaps.append((km, kmi))

def unregister():
    # Note: when unregistering, it's usually good practice to do it in reverse order you registered.
    # Can avoid strange issues like keymap still referring to operators already unregistered...
    # handle the keymap
    for km, kmi in addon_keymaps:
        km.keymap_items.remove(kmi)
    addon_keymaps.clear()

    bpy.utils.unregister_class(ObjectCursorArray)
    bpy.types.VIEW3D_MT_object.remove(menu_func)


if __name__ == "__main__":
    register()
