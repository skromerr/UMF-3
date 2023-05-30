using System;

namespace UMF_3;

public class FEM
{
   public delegate double Basis(Point3D point);

   private Grid grid;
   private Basis[] basis;
   private SparseMatrix globalMatrix = default!;
   private Vector globalVector = default!;
   private Vector solution = default!;
   private SLAE slae = default!;
   private Vector localVector, localVector1, localVector2;
   private Matrix stiffnessMatrix;
   private Matrix massMatrix;
   Test test = default!;

   public FEM(Grid Grid)
   {
      grid = Grid;
      stiffnessMatrix = new(8);
      massMatrix = new(8);
      localVector = new(16);
      localVector1 = new(8);
      localVector2 = new(8);
      basis = new Basis[] { ThreelinearBasis.Psi1, ThreelinearBasis.Psi2, ThreelinearBasis.Psi3, ThreelinearBasis.Psi4,
         ThreelinearBasis.Psi5, ThreelinearBasis.Psi6, ThreelinearBasis.Psi7, ThreelinearBasis.Psi8};
   }

   public void SetTest(Test test)
   {
      this.test = test;
   }

   public void SetSolver(SLAE slae)
   {
      this.slae = slae;
   }

   public void Compute()
   {
      BuildPortrait();

      AssemblySlae();
      AccountFirstConditions();

      slae.SetSLAE(globalVector, globalMatrix);
      solution = slae.Solve();

      PrintSolution();
   }

   public void BuildPortrait()
   {
      HashSet<int>[] list = new HashSet<int>[2*grid.Nodes.Length].Select(_ => new HashSet<int>()).ToArray();
      foreach (var element in grid.Elements)
      {
         foreach (var pos in element)
         {
            foreach (var node in element)
            {
               if (pos == node)
               {
                  list[2 * pos + 1].Add(2 * pos);
               }
               else if (pos > node)
               {
                  list[2 * pos].Add(2 * node);
                  list[2 * pos].Add(2 * node + 1);

                  list[2 * pos + 1].Add(2 * node);
                  list[2 * pos + 1].Add(2 * node + 1);
               }
            }
         }
      }

      list = list.Select(childlist => childlist.Order().ToHashSet()).ToArray();
      int count = list.Sum(childList => childList.Count);

      globalMatrix = new(2 * grid.Nodes.Length, count);
      globalVector = new(2 * grid.Nodes.Length);

      globalMatrix.Ig[0] = 0;

      for (int i = 0; i < list.Length; i++)
         globalMatrix.Ig[i + 1] = globalMatrix.Ig[i] + list[i].Count;

      int k = 0;

      foreach (var childList in list)
      {
         foreach (var value in childList)
         {
            globalMatrix.Jg[k++] = value;
         }
      }
   }

   public void AccountFirstConditions()
   {
      foreach (var node in grid.Boundaries) 
      {
         int row = 2 * node;

         globalMatrix.Di[row] = 1;
         globalVector[row] = test.Us(grid.Nodes[node]);

         for (int i = globalMatrix.Ig[row]; i < globalMatrix.Ig[row + 1]; i++)
         {
            globalMatrix.Ggl[i] = 0;
         }

         for (int col = row + 1; col < globalMatrix.Size; col++)
         {
            for (int j = globalMatrix.Ig[col]; j < globalMatrix.Ig[col + 1]; j++)
               if (globalMatrix.Jg[j] == row)
               {
                  globalMatrix.Ggu[j] = 0;
                  break;
               }
         }

         row = 2 * node + 1;

         globalMatrix.Di[row] = 1;
         globalVector[row] = test.Uc(grid.Nodes[node]);

         for (int i = globalMatrix.Ig[row]; i < globalMatrix.Ig[row + 1]; i++)
         {
            globalMatrix.Ggl[i] = 0;
         }

         for (int col = row + 1; col < globalMatrix.Size; col++)
         {
            for (int j = globalMatrix.Ig[col]; j < globalMatrix.Ig[col + 1]; j++)
               if (globalMatrix.Jg[j] == row)
               {
                  globalMatrix.Ggu[j] = 0;
                  break;
               }
         }
      }
   }

   private void AddElement(int i, int j, double value)
   {
      if (i == j)
      {
         globalMatrix.Di[i] += value;
         return;
      }

      if (i > j)
         for (int icol = globalMatrix.Ig[i]; icol < globalMatrix.Ig[i + 1]; icol++)
         {
            if (globalMatrix.Jg[icol] == j)
            {
               globalMatrix.Ggl[icol] += value;
               return;
            }
         }
      else
         for (int icol = globalMatrix.Ig[j]; icol < globalMatrix.Ig[j + 1]; icol++)
         {
            if (globalMatrix.Jg[icol] == i)
            {
               globalMatrix.Ggu[icol] += value;
               return;
            }
         }
   }

   private void AssemblySlae()
   {
      globalVector.Fill(0);
      globalMatrix.Clear();

      for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
      {
         AssemblyLocalSLAE(ielem);

         for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
            {
               AddElement(2 * grid.Elements[ielem][i], 2 * grid.Elements[ielem][j], stiffnessMatrix[i, j]);
               AddElement(2 * grid.Elements[ielem][i] + 1, 2 * grid.Elements[ielem][j] + 1, stiffnessMatrix[i, j]);
               AddElement(2 * grid.Elements[ielem][i], 2 * grid.Elements[ielem][j] + 1, -massMatrix[i, j]);
               AddElement(2 * grid.Elements[ielem][i] + 1, 2 * grid.Elements[ielem][j], massMatrix[i, j]);
            }

         AddElementToVector(ielem);

         stiffnessMatrix.Clear();
         massMatrix.Clear();
         localVector1.Fill(0);
         localVector2.Fill(0);
      }
   }

   private void AddElementToVector(int ielem)
   {
      for (int i = 0; i < 8; i++)
      {
         globalVector[2 * grid.Elements[ielem][i]] += localVector1[i];
         globalVector[2 * grid.Elements[ielem][i] + 1] += localVector2[i];
      }
   }

   private void AssemblyLocalSLAE(int ielem)
   {
      double hx = Math.Abs(grid.Nodes[grid.Elements[ielem][7]].X - grid.Nodes[grid.Elements[ielem][0]].X);
      double hy = Math.Abs(grid.Nodes[grid.Elements[ielem][7]].Y - grid.Nodes[grid.Elements[ielem][0]].Y);
      double hz = Math.Abs(grid.Nodes[grid.Elements[ielem][7]].Z - grid.Nodes[grid.Elements[ielem][0]].Z);

      double[,] matrixG = { { 1.0, -1.0 }, { -1.0, 1.0 } };
      double[,] matrixM = { { 2.0 / 6, 1.0 / 6 }, { 1.0 / 6, 2.0 / 6 } };

      for (int i = 0; i < 8; i++)
      {
         for (int j = 0; j < 8; j++)
         {
            stiffnessMatrix[i, j] = matrixG[GetMu(i), GetMu(j)] * matrixM[GetNu(i), GetNu(j)] * matrixM[GetTheta(i), GetTheta(j)] * hy * hz / hx +
               matrixM[GetMu(i), GetMu(j)] * matrixG[GetNu(i), GetNu(j)] * matrixM[GetTheta(i), GetTheta(j)] * hx * hz / hy +
               matrixM[GetMu(i), GetMu(j)] * matrixM[GetNu(i), GetNu(j)] * matrixG[GetTheta(i), GetTheta(j)] * hx * hy / hz;

            massMatrix[i, j] = matrixM[GetMu(i), GetMu(j)] * matrixM[GetNu(i), GetNu(j)] * matrixM[GetTheta(i), GetTheta(j)] * hx * hy * hz;
         }
      }

      for (int i = 0; i < 8; i++)
      {
         localVector1[i] = test.Fs(grid.Nodes[grid.Elements[ielem][i]]);
         localVector2[i] = test.Fc(grid.Nodes[grid.Elements[ielem][i]]);
      }

      localVector1 = massMatrix * localVector1;
      localVector2 = massMatrix * localVector2;

      stiffnessMatrix = grid.Lambda * stiffnessMatrix + (-grid.Omega) * grid.Omega * grid.Chi * massMatrix;
      massMatrix = grid.Omega * grid.Sigma * massMatrix;

      int GetMu(int i) => i % 2;
      int GetNu(int i) => i / 2 % 2;
      int GetTheta(int i) => i / 4;
   }

   public void PrintSolution()
   {
      Console.WriteLine(GetSolverName(slae.GetType().ToString()));
      //Console.WriteLine("Численное решение");
      for (int i = 0; i < solution.Length; i++)
      {
         //if (!grid.Boundaries.Contains(i / 2))
         //   Console.WriteLine(solution[i]);
      }
      Vector exactSolution = new(2 * grid.Nodes.Length);
      //Console.WriteLine("Точное решение");
      for (int i = 0; i < exactSolution.Length / 2; i++)
      {
         exactSolution[2 * i] = test.Us(grid.Nodes[i]);
         exactSolution[2 * i + 1] = test.Uc(grid.Nodes[i]);
         //if (!grid.Boundaries.Contains(i / 2))
         //   Console.WriteLine($"{exactSolution[2 * i]}\n{exactSolution[2 * i + 1]}");
      }
      //Console.WriteLine("Погрешность");
      Vector inaccuracySin = new(grid.Nodes.Length);
      Vector inaccuracyCos = new(grid.Nodes.Length);
      for (int i = 0; i < inaccuracySin.Length; i++)
      {
         if (!grid.Boundaries.Contains(i / 2))
         {
            inaccuracySin[i] = exactSolution[2 * i] - solution[2 * i];
            inaccuracyCos[i] = exactSolution[2 * i + 1] - solution[2 * i + 1];
            //Console.WriteLine($"{inaccuracy[i]}");
         }
      }
      Console.WriteLine("Относительная погрешность синус-компоненты");
      Console.WriteLine($"{inaccuracySin.Norm()/exactSolution.Norm()}");
      Console.WriteLine("Относительная погрешность косинус-компоненты");
      Console.WriteLine($"{inaccuracyCos.Norm() / exactSolution.Norm()}");
      Console.WriteLine($"Время: {slae.time}");
      Console.WriteLine($"Количество итераций : {slae.lastIter}");
      Console.WriteLine();
   }

   public void PrintSolution(StreamWriter sw)
   {
      sw.WriteLine(GetSolverName(slae.GetType().ToString()));
      //Console.WriteLine("Численное решение");
      for (int i = 0; i < solution.Length; i++)
      {
         //if (!grid.Boundaries.Contains(i / 2))
         //   Console.WriteLine(solution[i]);
      }
      Vector exactSolution = new(2 * grid.Nodes.Length);
      //Console.WriteLine("Точное решение");
      for (int i = 0; i < exactSolution.Length / 2; i++)
      {
         exactSolution[2 * i] = test.Us(grid.Nodes[i]);
         exactSolution[2 * i + 1] = test.Uc(grid.Nodes[i]);
         //if (!grid.Boundaries.Contains(i / 2))
         //   Console.WriteLine($"{exactSolution[2 * i]}\n{exactSolution[2 * i + 1]}");
      }
      //Console.WriteLine("Погрешность");
      Vector inaccuracySin = new(grid.Nodes.Length);
      Vector inaccuracyCos = new(grid.Nodes.Length);
      for (int i = 0; i < inaccuracySin.Length; i++)
      {
         if (!grid.Boundaries.Contains(i / 2))
         {
            inaccuracySin[i] = exactSolution[2 * i] - solution[2 * i];
            inaccuracyCos[i] = exactSolution[2 * i + 1] - solution[2 * i + 1];
            //Console.WriteLine($"{inaccuracy[i]}");
         }
      }
      sw.WriteLine("Относительная погрешность синус-компоненты");
      sw.WriteLine($"{inaccuracySin.Norm() / exactSolution.Norm()}");
      sw.WriteLine("Относительная погрешность косинус-компоненты");
      sw.WriteLine($"{inaccuracyCos.Norm() / exactSolution.Norm()}");
      sw.WriteLine($"Время: {slae.time}");
      sw.WriteLine($"Количество итераций : {slae.lastIter}");
      sw.WriteLine();
   }

   private string GetSolverName(string solverName)
   {
      string result = solverName[6..];
      result = result.Replace("Solver", ":");
      return result;
   }
}
