using UMF_3;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
Grid grid = new("C:/Users/Skromer/source/repos/UMF 3/UMF 3/grid.txt");
FEM test = new(grid);
test.SetTest(new Test1(grid));
test.SetSolver(new LUSolver());
test.Compute();

test.SetSolver(new LOSSolver());
test.Compute();

Console.WriteLine();