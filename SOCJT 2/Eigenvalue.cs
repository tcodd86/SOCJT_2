﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ConsoleApplication1
{
    class Eigenvalue
    {
        #region Properties
        private decimal npJ;
        public decimal pJ
        {
            get { return npJ; }
            set { npJ = value; }
        }

        private int nnJ;
        public int nJ
        {
            get { return nnJ; }
            set { nnJ = value; }
        }

        private decimal nSig;
        public decimal Sig
        {
            get { return nSig; }
            set { nSig = value; }
        }

        private double nEv;
        public double Ev
        {
            get { return nEv; }
            set { nEv = value; }
        }

        public bool isa1 { get; private set; }

        #endregion

        public Eigenvalue(decimal pJa, int nJa, decimal Siga, double Eva, bool isa1)
        {
            pJ = pJa;
            nJ = nJa;
            Sig = Siga;
            Ev = Eva;
            this.isa1 = isa1;
        }//end constructor

        /*
        public static Eigenvalue[] setAndSortEVs(List<double[]> evs, decimal S, bool inclSO)
        {
            List<Eigenvalue> eigen = new List<Eigenvalue>();
            int counter = 0;
            decimal J = 0.5M;
            S = S * -1M;
            decimal tempS = S;
            decimal maxS = S;
            if (inclSO == true)
            {
                maxS = maxS * -1M;
            }
            for (int i = 0; i < evs.Count; i++)
            {
                for (int j = 0; j < evs[i].Length; j++)
                {
                    eigen.Add(new Eigenvalue(J, j + 1, tempS, evs[i][j]));
                }
                if (tempS < maxS)
                {
                    tempS++;
                }
                else
                {
                    tempS = S;
                    J++;                    
                }
                counter += evs[i].Length;
            }
            Eigenvalue[] eigenarray = eigen.ToArray();
            bubbleSort(ref eigenarray);
            double ZPE = eigenarray[0].Ev;
            int[] temp = new int[evs.Count];
            for (int i = 0; i < evs.Count; i++)
            { 
                temp[i] = 1;
            }
            for (int i = 0; i < eigenarray.Length; i++)
            {
                eigenarray[i].Ev = eigenarray[i].Ev - ZPE;
            }
            int SOnumb = (int)(-2M * S) + 1;
            if (inclSO == false)
            {
                SOnumb = 1;
            }
            for (int i = 0; i < eigenarray.Length; i++)
            {
                int Snumb = (int)(eigenarray[i].Sig - S);
                int j = (int)(eigenarray[i].pJ - 0.5M);
                int place = j * SOnumb + Snumb;
                eigenarray[i].nJ = temp[place];
                temp[place]++;
            }
            return eigenarray;
        }

        private static void bubbleSort(ref Eigenvalue[] arr)
        {
            bool swapped = true;
            int j = 0;
            Eigenvalue tmp;
            while (swapped == true)
            {
                swapped = false;
                j++;
                for (int i = 0; i < arr.Length - j; i++)
                {
                    if (arr[i].Ev > arr[i + 1].Ev)
                    {
                        tmp = arr[i];
                        arr[i] = arr[i + 1];
                        arr[i + 1] = tmp;
                        swapped = true;
                    }//end if
                }//end for
            }//end while           
        }//end method bublleSort
        */

    }//end class Eigenvalue
}
