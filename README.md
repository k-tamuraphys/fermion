# fermion
For generation Hamiltonian with fermion operator

vm_dot(v,m)
  ベクトルvと行列mの内積を計算します. 結果はmと同じ型の行列です.
  
linear_tf(t,m)
  tを行列 m を各成分が行列であるような1次元配列(ベクトル)とします.
  t = (t_ij), m = (m_i)とした時,
  \sum_{j}t_ij m_j
  を返します.
  
mm_dot(m1, m2)
  m1 = (m1)_i m2 = (m2)_iを1次元配列とみなした時の内積を計算します.
  1つ1つの要素の積は行列同士の積です.
  
quad(m1, t, m2)
  m1 = (m1)_i, m2 = (m2)_i t=(t)_ijとした時に
  \sum_{ij}(m1)_i t_ij (m2)_j
  を返します.
  
  
fermionについて
  fermion.fermion(d = 1粒子状態の次元, 内部自由度の数)
    これで1粒子の状態空間を定義する.
  fermion.fermion(d, spin).c(i,alpha)
    i番目の, 内部自由度alphaでラベル付けされる状態のfermionの消滅演算子を生成
  fermion.fermion(d, spin).cdag(i,alpha)
    i番目の, 内部自由度alphaでラベル付けされる状態のfermionの生成演算子を生成
  fermion.fermion(d, spin).c_vec()
    [c(0,0), c(0,1), ... c(0,spin-1), c(1,0), c(1,1)...]の順に行列からなる1次元配列を生成
  fermion.fermion(d, spin).c_vec2()
    [c(0,0)c(0,0), c(0,0)c(0,1) ... c(0,0)c(dim-1, spin-1), c(0,1)c(0,0)...]の順に2粒子状態に対応した行列を1次元配列として生成
  fermion.fermion(d, spin).cdag_vec()
    c_vecのダガーを取ったもの
  fermion.fermion(d,spin).cdag_vec2()
    c_vec2のダガーを取ったもの
  fermion.fermion(d,spin).n_vec()
    c_vecと同じ順序でn(i,alpha)を並べた1次元配列を生成
  fermion.fermion(d, spin).h1(t_mat)
    t_matをhopping行列とするような1体のハミルトニアンを生成
  fermion.fermion(d, spin).h2(U_mat)
    U_matを相互作用の行列とするような2体のハミルトニアンを生成
