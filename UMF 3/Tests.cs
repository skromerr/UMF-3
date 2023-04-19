namespace UMF_3;

public abstract class Test
{
   protected double omega;
   protected double lambda;
   protected double sigma;
   protected double chi;

   public Test(Grid grid)
   {
      omega = grid.Omega;
      lambda = grid.Lambda;
      sigma = grid.Sigma;
      chi = grid.Chi;
   }

   public abstract double Us(Point3D point);
   protected abstract double divGradUs(Point3D point);
   public abstract double Uc(Point3D point);
   protected abstract double divGradUc(Point3D point);
   public double Fs(Point3D point) 
      => -lambda * divGradUs(point) - omega * sigma * Uc(point) - omega * omega * chi * Us(point);
   public double Fc(Point3D point)
      => -lambda * divGradUc(point) + omega * sigma * Us(point) - omega * omega * chi * Uc(point);
}

public class Test1 : Test
{
   public Test1(Grid grid) : base(grid) { }

   public override double Us(Point3D point) 
      => point.X + point.Y + point.Z;
   public override double Uc(Point3D point) 
      => point.Z - point.X - point.Y;
   protected override double divGradUs(Point3D point)
      => 0;
   protected override double divGradUc(Point3D point)
      => 0;
}