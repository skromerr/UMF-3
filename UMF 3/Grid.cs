using System.Drawing;

namespace UMF_3;

public class Point3D
{
   public double X { get; set; }
   public double Y { get; set; }
   public double Z { get; set; }

   public Point3D(double x, double y, double z)
   {
      X = x;
      Y = y;
      Z = z;
   }

   public static Point3D operator +(Point3D point, (double, double, double) value)
        => new(point.X + value.Item1, point.Y + value.Item2, point.Z + value.Item3);

   public static Point3D operator -(Point3D point, (double, double, double) value)
       => new(point.X - value.Item1, point.Y - value.Item2, point.Z - value.Item3);
}

public class Grid
{
   private double xStart;
   private double xEnd;
   private int xSteps;
   private double xRaz;
   private double yStart;
   private double yEnd;
   private int ySteps;
   private double yRaz;
   private double zStart;
   private double zEnd;
   private int zSteps;
   private double zRaz;
   public Point3D[] Nodes { get; private set; }
   public HashSet<int> Boundaries { get; private set; } 
   public int[][] Elements { get; private set; }
   public double Lambda { get; set; }
   public double Sigma { get; set; }
   public double Chi { get; set; }
   public double Omega { get; set; }

   public Grid(string path)
   {
      using (StreamReader sr = new (path))
      {
         string[] data;
         data = sr.ReadLine()!.Split(" ").ToArray();
         xStart = Convert.ToDouble(data[0]);
         xEnd = Convert.ToDouble(data[1]);
         xSteps = Convert.ToInt32(data[2]);
         xRaz = Convert.ToDouble(data[3]);

         data = sr.ReadLine()!.Split(" ").ToArray();
         yStart = Convert.ToDouble(data[0]);
         yEnd = Convert.ToDouble(data[1]);
         ySteps = Convert.ToInt32(data[2]);
         yRaz = Convert.ToDouble(data[3]);

         data = sr.ReadLine()!.Split(" ").ToArray();
         zStart = Convert.ToDouble(data[0]);
         zEnd = Convert.ToDouble(data[1]);
         zSteps = Convert.ToInt32(data[2]);
         zRaz = Convert.ToDouble(data[3]);

         data = sr.ReadLine()!.Split(" ").ToArray();
         Lambda = Convert.ToDouble(data[0]);
         Sigma = Convert.ToDouble(data[1]);
         Chi = Convert.ToDouble(data[2]);
         Omega = Convert.ToDouble(data[3]);
      }

      Elements = new int[xSteps * ySteps * zSteps].Select(_ => new int[8]).ToArray();
      Nodes = new Point3D[(xSteps + 1) * (ySteps + 1) * (zSteps + 1)];

      double sumRazX = 0, sumRazY = 0, sumRazZ = 0;
      for (int i = 0; i < xSteps; i++)
      {
         sumRazX += Math.Pow(xRaz, i);
      }

      for (int i = 0; i < ySteps; i++)
      {
         sumRazY += Math.Pow(yRaz, i);
      }

      for (int i = 0; i < zSteps; i++)
      {
         sumRazZ += Math.Pow(zRaz, i);
      }

      int nodesInRow = xSteps + 1;
      int nodesInSlice = nodesInRow * (ySteps + 1);

      double x = xStart, y = yStart, z = zStart;
      double xStep = (xEnd - xStart) / sumRazX, yStep = (yEnd - yStart) / sumRazY, zStep = (zEnd - zStart) / sumRazZ;

      Boundaries = new();

      for (int j = 0; j < xSteps; j++)
      {
         Nodes[j] = new(x, y, z);
         x += xStep;
         xStep *= xRaz;
         Boundaries.Add(j);
      }
      Nodes[xSteps] = new(xEnd, y, z);
      Boundaries.Add(xSteps);

      for (int i = 1; i < ySteps; i++)
      {
         y += yStep;
         yStep *= yRaz;
         for (int j = 0; j < xSteps + 1; j++)
         {
            Nodes[i * nodesInRow + j] = new(Nodes[j].X, y, z);
            Boundaries.Add(i * nodesInRow + j);
         }
      }

      for (int j = 0; j < xSteps + 1; j++)
      {
         Nodes[ySteps * nodesInRow + j] = new(Nodes[j].X, yEnd, z);
         Boundaries.Add(ySteps * nodesInRow + j);
      }

      for (int i = 1; i < zSteps; i++)
      {
         z += zStep;
         zStep *= zRaz;
         for (int j = 0; j < ySteps + 1; j++)
            for (int k = 0; k < xSteps + 1; k++)
            {
               Nodes[i * nodesInSlice + j * nodesInRow + k] = new(Nodes[k].X, Nodes[j * nodesInRow].Y, z);
               if (Nodes[k].X == xStart || Nodes[k].X == xEnd || Nodes[j * nodesInRow].Y == yStart || Nodes[j * nodesInRow].Y == yEnd)
                  Boundaries.Add(i * nodesInSlice + j * nodesInRow + k);
            }
      }

      for (int j = 0; j < ySteps + 1; j++)
         for (int k = 0; k < xSteps + 1; k++)
         {
            Nodes[zSteps * nodesInSlice + j * nodesInRow + k] = new(Nodes[k].X, Nodes[j * nodesInRow].Y, zEnd);
            Boundaries.Add(zSteps * nodesInSlice + j * nodesInRow + k);
         }

      int index = 0;

      for (int k = 0; k < zSteps; k++)
         for (int i = 0; i < ySteps; i++)
            for (int j = 0; j < xSteps; j++)
            {
               Elements[index][0] = j + nodesInRow * i + nodesInSlice * k;
               Elements[index][1] = (j + 1) + nodesInRow * i + nodesInSlice * k;
               Elements[index][2] = j + nodesInRow * (i + 1) + nodesInSlice * k;
               Elements[index][3] = (j + 1) + nodesInRow * (i + 1) + nodesInSlice * k;
               Elements[index][4] = j + nodesInRow * i + nodesInSlice * (k + 1);
               Elements[index][5] = (j + 1) + nodesInRow * i + nodesInSlice * (k + 1);
               Elements[index][6] = j + nodesInRow * (i + 1) + nodesInSlice * (k + 1);
               Elements[index++][7] = (j + 1) + nodesInRow * (i + 1) + nodesInSlice * (k + 1);
            }
   }
}
