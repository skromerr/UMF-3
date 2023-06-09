﻿using UMF_3;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");

double[] sigmas = { 10 };

Grid grid = new("C:/Users/Skromer/source/repos/UMF 3/UMF 3/grid.txt");

//StreamWriter sw = new("omega20k.txt", true);
//sw.Close();

foreach (var omega in sigmas)
{
   //sw = new("omega20k.txt", true);
   grid.Sigma = omega;
   Console.WriteLine($"Omega = {omega}");
   //sw.WriteLine($"Omega = {omega}");
   FEM test = new(grid);
   test.SetTest(new Test1(grid));
   //test.SetSolver(new LUSolver());
   //test.Compute();
   //test.PrintSolution(sw);

   test.SetSolver(new LOSSolver(1e-16, 2000));
   test.Compute();
   //test.PrintSolution(sw);
   Console.WriteLine("-------------------");
   //sw.WriteLine("---------------------");
   //sw.Close();
}

//sw.Close();

//SparseMatrix Matrix = new(5, 5);
//Matrix.Di = new double[] { 10, 20, 30, 40, 50 };
//Matrix.Ig = new int[] { 0, 0, 1, 2, 4, 5 };
//Matrix.Jg = new int[] { 0, 1, 0, 2, 2 };
//Matrix.Ggl = new double[] { 2, 3, 4, 8, 4 };
//Matrix.Ggu = new double[] { 5, 1, 2, 7, 3 };

//Vector vector = new(5);
//vector[0] = 28;
//vector[1] = 45;
//vector[2] = 139;
//vector[3] = 188;
//vector[4] = 262;

//SLAE slae = new LOSLUsqSolver();
//slae.SetSLAE(vector, Matrix);
//slae.Solve();
//slae.PrintSolution();

Console.WriteLine();