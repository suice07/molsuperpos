/*
 * @Author: Auier qi.mei@outlook.com
 * @Date: 2023-07-13 10:25:19
 * @LastEditors: Auier qi.mei@outlook.com
 * @LastEditTime: 2024-08-20 10:28:12
 * Copyright (c) 2024 by Auier qi.mei@outlook.com, All Rights Reserved. 
 */
#include <fstream>
#include <iostream>
#include <string>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/parsmart.h>
#include <openbabel/generic.h>

#include <Eigen/Dense>

using namespace std;
using namespace OpenBabel;
using namespace Eigen;

bool Superpose(OBMol& mol, OBMol& ref, string smarts);

int main(int argc, char* argv[]) {
   // setting forcefield
  string ff = "MMFF94";
  OBForceField* pFF = OBForceField::FindForceField(ff);
  if (!pFF) {
    cerr << "ForceField not found" << endl;
    return -1;
  }
  pFF->SetLogFile(&cout);
  pFF->SetLogLevel(OBFF_LOGLVL_NONE);

  // open input lig file to match
  string input_file = argv[1];
  ifstream ifs1(input_file);
  if (ifs1.fail()) {
    cerr << "input File not found" << endl;
    return -1;
  }
 // open reference lig file
  string ref_file = argv[2];
  ifstream ifs2(ref_file);
  if (ifs2.fail()) {
    cerr << "ref File not found" << endl;
    return -1;
  }

  // Loading reference molecules
  OBConversion conv;
  conv.SetInAndOutFormats("SDF", "SDF");
  OBMol mol, ref;
  conv.Read(&mol, &ifs1);
  conv.Read(&ref, &ifs2);
  ifs1.close();
  ifs2.close();

  // 两个分子的叠加计算
  string scaffold_smi = argv[3];
  Superpose(mol, ref, scaffold_smi);

  string out_file = argv[4];
  ofstream ofs(out_file);
  conv.Write(&mol, &ofs);
  ofs.close();

  return 0;
}

bool Superpose(OBMol& mol, OBMol& ref, string smarts) {
  // scaffold structure matching
  OBSmartsPattern smp;
  smp.Init(smarts);
  smp.Match(mol);
  vector<vector<int> > ml_mol = smp.GetUMapList();
  smp.Match(ref);
  vector<vector<int> > ml_ref = smp.GetUMapList();

  // Calculation of scaffold center coordinates
  Vector3d cm = Vector3d::Zero(3);
  Vector3d cr = Vector3d::Zero(3);
  for (int i = 0; i < ml_mol[0].size(); i++) {
    OBAtom* a_mol = mol.GetAtom(ml_mol[0][i]);
    cm(0) += a_mol->x();
    cm(1) += a_mol->y();
    cm(2) += a_mol->z();
    OBAtom* a_ref = ref.GetAtom(ml_ref[0][i]);
    cr(0) += a_ref->x();
    cr(1) += a_ref->y();
    cr(2) += a_ref->z();
  }
  cm /= (double)(ml_mol[0].size());
  cr /= (double)(ml_mol[0].size());

  // Move the scaffold center of two molecules to the origin
  for (auto ii = mol.BeginAtoms(); ii != mol.EndAtoms(); ii++) {
    double x = (*ii)->x() - cm(0);
    double y = (*ii)->y() - cm(1);
    double z = (*ii)->z() - cm(2);
    (*ii)->SetVector(x, y, z);
  }
  for (auto ii = ref.BeginAtoms(); ii != ref.EndAtoms(); ii++) {
    double x = (*ii)->x() - cr(0);
    double y = (*ii)->y() - cr(1);
    double z = (*ii)->z() - cr(2);
    (*ii)->SetVector(x, y, z);
  }

  // matrix calculation
  MatrixXd M = MatrixXd::Zero(4, 4);
  for (int i = 0; i < ml_mol[0].size(); i++) {
    OBAtom* a_mol = mol.GetAtom(ml_mol[0][i]);
    double x = a_mol->x();
    double y = a_mol->y();
    double z = a_mol->z();
    OBAtom* a_ref = ref.GetAtom(ml_ref[0][i]);
    double x0 = a_ref->x();
    double y0 = a_ref->y();
    double z0 = a_ref->z();
    M(0, 0) += -x * x0 - y * y0 - z * z0;
    M(1, 0) += z * y0 - y * z0;
    M(1, 1) += -x * x0 + y * y0 + z * z0;
    M(2, 0) += x * z0 - z * x0;
    M(2, 1) += -x * y0 - y * x0;
    M(2, 2) += x * x0 - y * y0 + z * z0;
    M(3, 0) += y * x0 - x * y0;
    M(3, 1) += -x * z0 - z * x0;
    M(3, 2) += -y * z0 - z * y0;
    M(3, 3) += x * x0 + y * y0 - z * z0;
  }

  // matrix diagonalization
  SelfAdjointEigenSolver<MatrixXd> es(M);
  if (es.info() != Success) exit(0);
  VectorXd E;
  MatrixXd V;
  E = es.eigenvalues();
  V = es.eigenvectors();

  // compute rotation matrix from quaternion
  Matrix3d R;
  double q0 = V(0, 0);
  double q1 = V(1, 0);
  double q2 = V(2, 0);
  double q3 = V(3, 0);
  R(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
  R(1, 0) = 2 * (q0 * q3 + q1 * q2);
  R(2, 0) = 2 * (q1 * q3 - q0 * q2);
  R(0, 1) = 2 * (q1 * q2 - q0 * q3);
  R(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
  R(2, 1) = 2 * (q2 * q3 + q0 * q1);
  R(0, 2) = 2 * (q0 * q2 + q1 * q3);
  R(1, 2) = 2 * (-q0 * q1 + q2 * q3);
  R(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

  // molecular rotation and translation
  for (auto ii = mol.BeginAtoms(); ii != mol.EndAtoms(); ii++) {
    Vector3d X;
    X(0) = (*ii)->x();
    X(1) = (*ii)->y();
    X(2) = (*ii)->z();
    X = R * X;   // rotate
    X = X + cr;  // translation
    (*ii)->SetVector(X(0), X(1), X(2));
  }

  return true;
}
