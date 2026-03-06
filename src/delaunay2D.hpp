#pragma once
// THE LIGHT VERSION FROM NR3

typedef unsigned long int ulint;
// for windows which is a 32 bits OS
// it may be nessesary to replace the line above by:
// typedef unsigned int ulint;

struct Point2D {
  double x, y;
  Point2D() : x(0.0), y(0.0) {}
  Point2D(double X, double Y) : x(X), y(Y) {}
};

double dist(const Point2D &p, const Point2D &q) {
  double dd = SQR(q.x - p.x) + SQR(q.y - p.y);
  return sqrt(dd);
}

struct Ranhash {
  ulint int64(ulint u) {
    ulint v = u * 3935559000370003845LL + 2691343689449507681LL;
    v ^= v >> 21;
    v ^= v << 37;
    v ^= v >> 4;
    v *= 4768777513237032717LL;
    v ^= v << 20;
    v ^= v >> 41;
    v ^= v << 5;
    return v;
  }
  inline unsigned int int32(ulint u) { return (unsigned int)(int64(u) & 0xffffffff); }
  inline double doub(ulint u) { return 5.42101086242752217E-20 * int64(u); }
};

// HASH -----

template <class keyT, class hfnT> struct Hashtable {
  hfnT hash;
  int nhash, nmax, nn, ng;
  vector<int> htable, next, garbg;
  vector<ulint> thehash;
  Hashtable(int nh, int nv);
  int iget(const keyT &key);
  int iset(const keyT &key);
  int ierase(const keyT &key);
  int ireserve();
  int irelinquish(int k);
};

template <class keyT, class hfnT>
Hashtable<keyT, hfnT>::Hashtable(int nh, int nv)
    : hash(sizeof(keyT)), nhash(nh), nmax(nv), nn(0), ng(0), htable(nh), next(nv), garbg(nv), thehash(nv) {
  for (int j = 0; j < nh; j++) {
    htable[j] = -1;
  }
}

template <class keyT, class hfnT> int Hashtable<keyT, hfnT>::iget(const keyT &key) {
  int j, k;
  ulint pp = hash.fn(&key);
  j = (int)(pp % nhash);
  for (k = htable[j]; k != -1; k = next[k]) {
    if (thehash[k] == pp) {
      return k;
    }
  }
  return -1;
}

template <class keyT, class hfnT> int Hashtable<keyT, hfnT>::iset(const keyT &key) {
  int j, k, kprev;
  ulint pp = hash.fn(&key);
  j = (int)(pp % nhash);
  if (htable[j] == -1) {
    k = ng ? garbg[--ng] : nn++;
    htable[j] = k;
  } else {
    for (k = htable[j]; k != -1; k = next[k]) {
      if (thehash[k] == pp) {
        return k;
      }
      kprev = k;
    }
    k = ng ? garbg[--ng] : nn++;
    next[kprev] = k;
  }
  if (k >= nmax)
    cerr << "storing too many values" << endl;
  thehash[k] = pp;
  next[k] = -1;
  return k;
}

template <class keyT, class hfnT> int Hashtable<keyT, hfnT>::ierase(const keyT &key) {
  int j, k, kprev;
  ulint pp = hash.fn(&key);
  j = (int)(pp % nhash);
  if (htable[j] == -1)
    return -1;
  kprev = -1;
  for (k = htable[j]; k != -1; k = next[k]) {
    if (thehash[k] == pp) {
      if (kprev == -1)
        htable[j] = next[k];
      else
        next[kprev] = next[k];
      garbg[ng++] = k;
      return k;
    }
    kprev = k;
  }
  return -1;
}

template <class keyT, class hfnT> int Hashtable<keyT, hfnT>::ireserve() {
  int k = ng ? garbg[--ng] : nn++;
  if (k >= nmax)
    cerr << "reserving too many values" << endl;
  next[k] = -2;
  return k;
}

template <class keyT, class hfnT> int Hashtable<keyT, hfnT>::irelinquish(int k) {
  if (next[k] != -2) {
    return -1;
  }
  garbg[ng++] = k;
  return k;
}

template <class keyT, class elT, class hfnT> struct Hash : Hashtable<keyT, hfnT> {
  using Hashtable<keyT, hfnT>::iget;
  using Hashtable<keyT, hfnT>::iset;
  using Hashtable<keyT, hfnT>::ierase;
  vector<elT> els;

  Hash(int nh, int nm) : Hashtable<keyT, hfnT>(nh, nm), els(nm) {}

  void set(const keyT &key, const elT &el) { els[iset(key)] = el; }

  int get(const keyT &key, elT &el) {
    int ll = iget(key);
    if (ll < 0)
      return 0;
    el = els[ll];
    return 1;
  }

  elT &operator[](const keyT &key) {
    int ll = iget(key);
    if (ll < 0) {
      ll = iset(key);
      els[ll] = elT();
    }
    return els[ll];
  }

  int count(const keyT &key) {
    int ll = iget(key);
    return (ll < 0 ? 0 : 1);
  }

  int erase(const keyT &key) { return (ierase(key) < 0 ? 0 : 1); }
};

//  ----------

struct Circle {
  Point2D center;
  double radius;
  Circle(const Point2D &cen, double rad) : center(cen), radius(rad) {}
};

Circle circumcircle(Point2D a, Point2D b, Point2D c) {
  double a0, a1, c0, c1, det, asq, csq, ctr0, ctr1, rad2;
  a0 = a.x - b.x;
  a1 = a.y - b.y;
  c0 = c.x - b.x;
  c1 = c.y - b.y;
  det = a0 * c1 - c0 * a1;
  if (det == 0.0)
    cerr << "no circle thru colinear points" << endl;
  det = 0.5 / det;
  asq = a0 * a0 + a1 * a1;
  csq = c0 * c0 + c1 * c1;
  ctr0 = det * (asq * c1 - csq * a1);
  ctr1 = det * (csq * a0 - asq * c0);
  rad2 = ctr0 * ctr0 + ctr1 * ctr1;
  return Circle(Point2D(ctr0 + b.x, ctr1 + b.y), sqrt(rad2));
}

//  ----------

struct Triel {
  Point2D *pts;
  int p[3];
  int d[3];
  int stat;
  void setme(int a, int b, int c, Point2D *ptss) {
    pts = ptss;
    p[0] = a;
    p[1] = b;
    p[2] = c;
    d[0] = d[1] = d[2] = -1;
    stat = 1;
  }
  int contains(Point2D point) {
    double d;
    int i, j, ztest = 0;
    for (i = 0; i < 3; i++) {
      j = (i + 1) % 3;
      d = (pts[p[j]].x - pts[p[i]].x) * (point.y - pts[p[i]].y) -
          (pts[p[j]].y - pts[p[i]].y) * (point.x - pts[p[i]].x);
      if (d < 0.0)
        return -1;
      if (d == 0.0)
        ztest = 1;
    }
    return (ztest ? 0 : 1);
  }
};
double incircle(Point2D d, Point2D a, Point2D b, Point2D c) {
  Circle cc = circumcircle(a, b, c);
  double radd = SQR(d.x - cc.center.x) + SQR(d.y - cc.center.y);
  return (SQR(cc.radius) - radd);
}

struct Nullhash {
  Nullhash(int) {}
  ulint fn(const void *key) const { return *((ulint *)key); }
};

struct Delaunay {
  int npts, ntri, ntree, ntreemax, opt;
  double delx, dely;
  vector<Point2D> pts;
  vector<Triel> thelist;
  Hash<ulint, int, Nullhash> *linehash;
  Hash<ulint, int, Nullhash> *trihash;
  int *perm;
  Delaunay(vector<Point2D> &pvec, int options = 0);
  Ranhash hashfn;
  double interpolate(const Point2D &p, const vector<double> &fnvals, double defaultval = 0.0);
  void insertapoint(int r);
  int whichcontainspt(const Point2D &p, int strict = 0);
  int storetriangle(int a, int b, int c);
  void erasetriangle(int a, int b, int c, int d0, int d1, int d2);
  static unsigned int jran;
  static const double fuzz, bigscale;
};

const double Delaunay::fuzz = 1.0e-6;
const double Delaunay::bigscale = 1000.0;
unsigned int Delaunay::jran = 14921620;

Delaunay::Delaunay(vector<Point2D> &pvec, int options)
    : npts(pvec.size()), ntri(0), ntree(0), ntreemax(10 * npts + 1000), opt(options), pts(npts + 3),
      thelist(ntreemax) {
  int j;
  double xl, xh, yl, yh;
  linehash = new Hash<ulint, int, Nullhash>(6 * npts + 12, 6 * npts + 12);
  trihash = new Hash<ulint, int, Nullhash>(2 * npts + 6, 2 * npts + 6);
  perm = new int[npts];
  xl = xh = pvec[0].x;
  yl = yh = pvec[0].y;
  for (j = 0; j < npts; j++) {
    pts[j] = pvec[j];
    perm[j] = j;
    if (pvec[j].x < xl)
      xl = pvec[j].x;
    if (pvec[j].x > xh)
      xh = pvec[j].x;
    if (pvec[j].y < yl)
      yl = pvec[j].y;
    if (pvec[j].y > yh)
      yh = pvec[j].y;
  }
  delx = xh - xl;
  dely = yh - yl;
  pts[npts] = Point2D(0.5 * (xl + xh), yh + bigscale * dely);
  pts[npts + 1] = Point2D(xl - 0.5 * bigscale * delx, yl - 0.5 * bigscale * dely);
  pts[npts + 2] = Point2D(xh + 0.5 * bigscale * delx, yl - 0.5 * bigscale * dely);
  storetriangle(npts, npts + 1, npts + 2);
  for (j = npts; j > 0; j--)
    SWAP(perm[j - 1], perm[hashfn.int64(jran++) % j]);
  for (j = 0; j < npts; j++)
    insertapoint(perm[j]);
  for (j = 0; j < ntree; j++) {
    if (thelist[j].stat > 0) {
      if (thelist[j].p[0] >= npts || thelist[j].p[1] >= npts || thelist[j].p[2] >= npts) {
        thelist[j].stat = -1;
        ntri--;
      }
    }
  }
  if (!(opt & 1)) {
    delete[] perm;
    delete trihash;
    delete linehash;
  }
}
void Delaunay::insertapoint(int r) {
  int i, j, k, l, s, tno, ntask, d0, d1, d2;
  ulint key;
  int tasks[50], taski[50], taskj[50];
  for (j = 0; j < 3; j++) {
    tno = whichcontainspt(pts[r], 1);
    if (tno >= 0)
      break;
    pts[r].x += fuzz * delx * (hashfn.doub(jran++) - 0.5);
    pts[r].y += fuzz * dely * (hashfn.doub(jran++) - 0.5);
  }
  if (j == 3)
    cerr << "points degenerate even after fuzzing" << endl;
  ntask = 0;
  i = thelist[tno].p[0];
  j = thelist[tno].p[1];
  k = thelist[tno].p[2];
  if (opt & 2 && i < npts && j < npts && k < npts)
    return;
  d0 = storetriangle(r, i, j);
  tasks[++ntask] = r;
  taski[ntask] = i;
  taskj[ntask] = j;
  d1 = storetriangle(r, j, k);
  tasks[++ntask] = r;
  taski[ntask] = j;
  taskj[ntask] = k;
  d2 = storetriangle(r, k, i);
  tasks[++ntask] = r;
  taski[ntask] = k;
  taskj[ntask] = i;
  erasetriangle(i, j, k, d0, d1, d2);
  while (ntask) {
    s = tasks[ntask];
    i = taski[ntask];
    j = taskj[ntask--];
    key = hashfn.int64(j) - hashfn.int64(i);
    if (!linehash->get(key, l))
      continue;
    if (incircle(pts[l], pts[j], pts[s], pts[i]) > 0.0) {
      d0 = storetriangle(s, l, j);
      d1 = storetriangle(s, i, l);
      erasetriangle(s, i, j, d0, d1, -1);
      erasetriangle(l, j, i, d0, d1, -1);
      key = hashfn.int64(i) - hashfn.int64(j);
      linehash->erase(key);
      key = 0 - key;
      linehash->erase(key);
      tasks[++ntask] = s;
      taski[ntask] = l;
      taskj[ntask] = j;
      tasks[++ntask] = s;
      taski[ntask] = i;
      taskj[ntask] = l;
    }
  }
}
int Delaunay::whichcontainspt(const Point2D &p, int strict) {
  int i, j, k = 0;
  while (thelist[k].stat <= 0) {
    for (i = 0; i < 3; i++) {
      if ((j = thelist[k].d[i]) < 0)
        continue;
      if (strict) {
        if (thelist[j].contains(p) > 0)
          break;
      } else {
        if (thelist[j].contains(p) >= 0)
          break;
      }
    }
    if (i == 3)
      return -1;
    k = j;
  }
  return k;
}

void Delaunay::erasetriangle(int a, int b, int c, int d0, int d1, int d2) {
  ulint key;
  int j = 0;
  key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
  if (trihash->get(key, j) == 0)
    cerr << "nonexistent triangle" << endl;
  trihash->erase(key);
  thelist[j].d[0] = d0;
  thelist[j].d[1] = d1;
  thelist[j].d[2] = d2;
  thelist[j].stat = 0;
  ntri--;
}

int Delaunay::storetriangle(int a, int b, int c) {
  ulint key;
  thelist[ntree].setme(a, b, c, &pts[0]);
  key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
  trihash->set(key, ntree);
  key = hashfn.int64(b) - hashfn.int64(c);
  linehash->set(key, a);
  key = hashfn.int64(c) - hashfn.int64(a);
  linehash->set(key, b);
  key = hashfn.int64(a) - hashfn.int64(b);
  linehash->set(key, c);
  if (++ntree == ntreemax)
    cerr << "thelist is sized too small" << endl;
  ntri++;
  return (ntree - 1);
}

double Delaunay::interpolate(const Point2D &p, const vector<double> &fnvals, double defaultval) {
  int n, i, j, k;
  double wgts[3];
  int ipts[3];
  double sum, ans = 0.0;
  n = whichcontainspt(p);
  if (n < 0)
    return defaultval;
  for (i = 0; i < 3; i++)
    ipts[i] = thelist[n].p[i];
  for (i = 0, j = 1, k = 2; i < 3; i++, j++, k++) {
    if (j == 3)
      j = 0;
    if (k == 3)
      k = 0;
    wgts[k] = (pts[ipts[j]].x - pts[ipts[i]].x) * (p.y - pts[ipts[i]].y) -
              (pts[ipts[j]].y - pts[ipts[i]].y) * (p.x - pts[ipts[i]].x);
  }
  sum = wgts[0] + wgts[1] + wgts[2];
  if (sum == 0)
    cerr << "degenerate triangle" << endl;
  for (i = 0; i < 3; i++)
    ans += wgts[i] * fnvals[ipts[i]] / sum;
  return ans;
}
