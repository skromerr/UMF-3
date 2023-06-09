﻿namespace UMF_3;
public static class ThreelinearBasis
{
   public static double Psi1(Point3D point)
      => (1 - point.X) * (1 - point.Y) * (1 - point.Z);
   public static double Psi2(Point3D point)
      => point.X * (1 - point.Y) * (1 - point.Z);
   public static double Psi3(Point3D point)
      => (1 - point.X) * point.Y * (1 - point.Z);
   public static double Psi4(Point3D point)
      => point.X * point.Y * (1 - point.Z);
   public static double Psi5(Point3D point)
      => (1 - point.X) * (1 - point.Y) * point.Z;
   public static double Psi6(Point3D point)
      => point.X * (1 - point.Y) * point.Z;
   public static double Psi7(Point3D point)
      => (1 - point.X) * point.Y * point.Z;
   public static double Psi8(Point3D point)
      => point.X * point.Y * point.Z;
}
