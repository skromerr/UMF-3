﻿using System.ComponentModel.Design;
using System.Diagnostics;
using System.Numerics;

namespace UMF_3;

public abstract class SLAE
{
   protected SparseMatrix matrix = default!;
   protected Vector vector = default!;
   protected Vector solution = default!;
   public double time;

   public void SetSLAE(Vector vector, SparseMatrix matrix)
   {
      this.vector = vector;
      this.matrix = matrix;
   }

   public abstract Vector Solve();


   protected void LU()
   {
      for (int i = 0; i < matrix.Size; i++)
      {

         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
         {
            int jCol = matrix.Jg[j];
            int jk = matrix.Ig[jCol];
            int k = matrix.Ig[i];

            int sdvig = matrix.Jg[matrix.Ig[i]] - matrix.Jg[matrix.Ig[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < matrix.Ig[jCol + 1]; k++, jk++)
            {
               sumL += matrix.Ggl[k] * matrix.Ggu[jk];
               sumU += matrix.Ggu[k] * matrix.Ggl[jk];
            }

            matrix.Ggl[j] -= sumL;
            matrix.Ggu[j] -= sumU;
            matrix.Ggu[j] /= matrix.Di[jCol];
         }

         double sumD = 0.0;
         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
            sumD += matrix.Ggl[j] * matrix.Ggu[j];

         matrix.Di[i] -= sumD;
      }
   }

   protected void ForwardElimination()
   {
      for (int i = 0; i < matrix.Size; i++)
      {
         for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
         {
            solution[i] -= matrix.Ggl[j] * solution[matrix.Jg[j]];
         }

         solution[i] /= matrix.Di[i];
      }
   }

   protected void BackwardSubstitution()
   {
      for (int i = matrix.Size - 1; i >= 0; i--)
      {
         for (int j = matrix.Ig[i + 1] - 1; j >= matrix.Ig[i]; j--)
         {
            solution[matrix.Jg[j]] -= matrix.Ggu[j] * solution[i];
         }
      }
   }

   public void PrintSolution()
   {
      for(int i = 0; i < solution.Length; i++)
      {
         Console.WriteLine(solution[i]);
      }
   }
}

public class LOSSolver : SLAE
{
   private double eps = 1e-16;
   private int maxIters = 2000;
   public override Vector Solve()
   {
      solution = new(vector.Length);
      Vector.Copy(vector, solution);


      Vector r = vector - matrix * solution;
      Vector z = 1 * r;
      Vector p = matrix * z;
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double nevyaz = r * r;


      for (int i = 1; i <= maxIters & nevyaz > eps; i++)
      {
         alpha = (p * r) / (p * p);
         solution += alpha * z;

         r -= alpha * p;
         tmp = matrix * r;

         beta = -(p * tmp) / (p * p);

         z = r + beta * z;
         p = tmp + beta * p;

         nevyaz = r * r;
      }

      sw.Stop();
      time = sw.ElapsedMilliseconds;

      return solution;
   }
}

public class LUSolver : SLAE
{
   public override Vector Solve()
   {
      solution = new(vector.Length);
      Vector.Copy(vector, solution);
      matrix = matrix.ConvertToProfile();

      Stopwatch sw = Stopwatch.StartNew();

      LU();
      ForwardElimination();
      BackwardSubstitution();

      sw.Stop();
      time = sw.ElapsedMilliseconds;

      return solution;
   }
}

public class LOSLUSolver : SLAE
{
   private double eps = 1e-16;
   private int maxIters = 2000;
   public override Vector Solve() 
   {
      solution = new(vector.Length);
      Vector.Copy(vector, solution);

      SparseMatrix matrixLU = new(matrix.Size, matrix.Jg.Length);
      SparseMatrix.Copy(matrix, matrixLU);

      //matrixLU = matrixLU.ConvertToProfile();

      LU(matrixLU);

      Vector r = DirElim(matrixLU, vector - matrix * solution);
      Vector z = BackSub(matrixLU, r);
      Vector p = DirElim(matrixLU, matrix * z);
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double nevyaz = r * r;


      for (int i = 1; i <= maxIters & nevyaz > eps; i++)
      {
         alpha = (p * r) / (p * p);
         solution += alpha * z;

         r -= alpha * p;

         tmp = DirElim(matrixLU, matrix * BackSub(matrixLU, r));
         beta = -(p * tmp) / (p * p);

         z = BackSub(matrixLU, r) + beta * z;
         p = tmp + beta * p;

         nevyaz = r * r;
      }

      sw.Stop();
      time = sw.ElapsedMilliseconds;

      return solution;
   }

   protected static void LU(SparseMatrix Matrix)
   {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
         {
            int jCol = Matrix.Jg[j];
            int jk = Matrix.Ig[jCol];
            int k = Matrix.Ig[i];

            int sdvig = Matrix.Jg[Matrix.Ig[i]] - Matrix.Jg[Matrix.Ig[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix.Ig[jCol + 1]; k++, jk++)
            {
               sumL += Matrix.Ggl[k] * Matrix.Ggu[jk];
               sumU += Matrix.Ggu[k] * Matrix.Ggl[jk];
            }

            Matrix.Ggl[j] -= sumL;
            Matrix.Ggu[j] -= sumU;
            Matrix.Ggu[j] /= Matrix.Di[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
            sumD += Matrix.Ggl[j] * Matrix.Ggu[j];

         Matrix.Di[i] -= sumD;
      }
   }

   protected static Vector DirElim(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Length);
      Vector.Copy(b, result);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
         {
            result[i] -= Matrix.Ggl[j] * result[Matrix.Jg[j]];
         }

         result[i] /= Matrix.Di[i];
      }

      return result;
   }

   protected static Vector BackSub(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Length);
      Vector.Copy(b, result);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         for (int j = Matrix.Ig[i + 1] - 1; j >= Matrix.Ig[i]; j--)
         {
            result[Matrix.Jg[j]] -= Matrix.Ggu[j] * result[i];
         }
      }

      return result;
   }

}

public class LOSLUsqSolver : SLAE
{
   protected double eps = 1e-16;
   protected int maxIters = 2000;
   public override Vector Solve() 
   {
      solution = new(vector.Length);
      Vector.Copy(vector, solution);

      SparseMatrix matrixLU = new(matrix.Size, matrix.Jg.Length);
      SparseMatrix.Copy(matrix, matrixLU);

      //matrixLU = matrixLU.ConvertToProfile();

      LUsq(matrixLU);

      Vector r = DirElim(matrixLU, vector - matrix * solution);
      Vector z = BackSub(matrixLU, r);
      Vector p = DirElim(matrixLU, matrix * z);
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double nevyaz = r * r;


      for (int i = 1; i <= maxIters & nevyaz > eps; i++)
      {
         alpha = (p * r) / (p * p);
         solution += alpha * z;

         r -= alpha * p;

         tmp = DirElim(matrixLU, matrix * BackSub(matrixLU, r));
         beta = -(p * tmp) / (p * p);

         z = BackSub(matrixLU, r) + beta * z;
         p = tmp + beta * p;

         nevyaz = r * r;
      }

      sw.Stop();
      time = sw.ElapsedMilliseconds;

      return solution;
   }

   protected static void LUsq(SparseMatrix Matrix)
   {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
         {
            int jCol = Matrix.Jg[j];
            int jk = Matrix.Ig[jCol];
            int k = Matrix.Ig[i];

            int sdvig = Matrix.Jg[Matrix.Ig[i]] - Matrix.Jg[Matrix.Ig[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix.Ig[jCol + 1]; k++, jk++)
            {
               sumL += Matrix.Ggl[k] * Matrix.Ggu[jk];
               sumU += Matrix.Ggu[k] * Matrix.Ggl[jk];
            }

            Matrix.Ggl[j] -= sumL;
            Matrix.Ggl[j] /= Matrix.Di[jCol];
            Matrix.Ggu[j] -= sumU;
            Matrix.Ggu[j] /= Matrix.Di[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
            sumD += Matrix.Ggl[j] * Matrix.Ggu[j];

         Matrix.Di[i] = Math.Sqrt(Matrix.Di[i] - sumD);
      }
   }

   protected static Vector DirElim(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Length);
      Vector.Copy(b, result);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
         {
            result[i] -= Matrix.Ggl[j] * result[Matrix.Jg[j]];
         }

         result[i] /= Matrix.Di[i];
      }

      return result;
   }

   protected static Vector BackSub(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Length);
      Vector.Copy(b, result);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         result[i] /= Matrix.Di[i];

         for (int j = Matrix.Ig[i + 1] - 1; j >= Matrix.Ig[i]; j--)
         {
            result[Matrix.Jg[j]] -= Matrix.Ggu[j] * result[i];
         }


      }

      return result;
   }

}

public class BCGSTABSolver : SLAE
{
   private double eps = 1e-16;
   private int maxIters = 2000;
   public override Vector Solve()
   {
      solution = new(vector.Length);
      Vector.Copy(vector, solution);

      SparseMatrix matrixLU = new(matrix.Size, matrix.Jg.Length);
      SparseMatrix.Copy(matrix, matrixLU);

      LU(matrixLU);


      
      return solution;
   }

   protected static void LU(SparseMatrix Matrix)
   {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
         {
            int jCol = Matrix.Jg[j];
            int jk = Matrix.Ig[jCol];
            int k = Matrix.Ig[i];

            int sdvig = Matrix.Jg[Matrix.Ig[i]] - Matrix.Jg[Matrix.Ig[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix.Ig[jCol + 1]; k++, jk++)
            {
               sumL += Matrix.Ggl[k] * Matrix.Ggu[jk];
               sumU += Matrix.Ggu[k] * Matrix.Ggl[jk];
            }

            Matrix.Ggl[j] -= sumL;
            Matrix.Ggu[j] -= sumU;
            Matrix.Ggu[j] /= Matrix.Di[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
            sumD += Matrix.Ggl[j] * Matrix.Ggu[j];

         Matrix.Di[i] -= sumD;
      }
   }

   protected static Vector DirElim(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Length);
      Vector.Copy(b, result);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix.Ig[i]; j < Matrix.Ig[i + 1]; j++)
         {
            result[i] -= Matrix.Ggl[j] * result[Matrix.Jg[j]];
         }

         result[i] /= Matrix.Di[i];
      }

      return result;
   }

   protected static Vector BackSub(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Length);
      Vector.Copy(b, result);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         for (int j = Matrix.Ig[i + 1] - 1; j >= Matrix.Ig[i]; j--)
         {
            result[Matrix.Jg[j]] -= Matrix.Ggu[j] * result[i];
         }
      }
      return result;
   }

}