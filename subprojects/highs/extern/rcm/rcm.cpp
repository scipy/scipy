#include "rcm.h"

#include <cstdio>
#include <cstdlib>

#include "OrderingPrint.h"

// offset from C to Fortran numbering
const HighsInt offset = 1;

//****************************************************************************80

static void Highs_rcm_degree(HighsInt root, HighsInt adj_num,
                             const HighsInt adj_row[], const HighsInt adj[],
                             HighsInt mask[], HighsInt deg[], HighsInt* iccsze,
                             HighsInt ls[], HighsInt node_num)

//****************************************************************************80
//
//  Purpose:
//
//    degree() computes the degrees of the nodes in the connected component.
//
//  Discussion:
//
//    The connected component is specified by MASK and ROOT.
//    Nodes for which MASK is zero are ignored.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    This version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Input:
//
//    int ROOT, the node that defines the connected component.
//
//    int ADJ_NUM, the number of adjacency entries.
//
//    int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    int MASK[NODE_NUM], is nonzero for those nodes which are
//    to be considered.
//
//  Output:
//
//    int DEG[NODE_NUM], contains, for each  node in the connected
//    component, its degree.
//
//    int *ICCSIZE, the number of nodes in the connected component.
//
//    int LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
//    in the connected component, starting with ROOT, and proceeding
//    by levels.
//
//    int NODE_NUM, the number of nodes.
//
{
  HighsInt i;
  HighsInt ideg;
  HighsInt j;
  HighsInt jstop;
  HighsInt jstrt;
  HighsInt lbegin;
  HighsInt lvlend;
  HighsInt lvsize;
  HighsInt nbr;
  HighsInt node;
  //
  //  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
  //  Use a separate array instead, to keep adj_row const.
  HighsInt* mark = new HighsInt[node_num];
  for (HighsInt i = 0; i < node_num; ++i) mark[i] = 1;

  ls[0] = root;
  mark[root - 1] = -1;
  lvlend = 0;
  *iccsze = 1;
  //
  //  LBEGIN is the pointer to the beginning of the current level, and
  //  LVLEND points to the end of this level.
  //
  for (;;) {
    lbegin = lvlend + 1;
    lvlend = *iccsze;
    //
    //  Find the degrees of nodes in the current level,
    //  and at the same time, generate the next level.
    //
    for (i = lbegin; i <= lvlend; i++) {
      node = ls[i - 1];
      jstrt = adj_row[node - 1] + offset;
      jstop = adj_row[node] - 1 + offset;
      ideg = 0;

      for (j = jstrt; j <= jstop; j++) {
        nbr = adj[j - 1] + offset;

        if (mask[nbr - 1] != 0) {
          ideg = ideg + 1;

          if (0 <= mark[nbr - 1]) {
            mark[nbr - 1] = -1;
            *iccsze = *iccsze + 1;
            ls[*iccsze - 1] = nbr;
          }
        }
      }
      deg[node - 1] = ideg;
    }
    //
    //  Compute the current level width.
    //
    lvsize = *iccsze - lvlend;
    //
    //  If the current level width is nonzero, generate another level.
    //
    if (lvsize == 0) {
      break;
    }
  }
  //
  //  Reset ADJ_ROW to its correct sign and return.
  //
  for (i = 0; i < *iccsze; i++) {
    node = ls[i] - 1;
  }

  delete[] mark;

  return;
}
//****************************************************************************80

static void Highs_rcm_i4vec_reverse(HighsInt n, HighsInt a[])

//****************************************************************************80
//
//  Purpose:
//
//    i4vec_reverse() reverses the elements of an I4VEC.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Input:
//
//    int N, the number of entries in the array.
//
//    int A(N), the array to be reversed.
//
//  Output:
//
//    int A[N]: the reversed array.
//
{
  HighsInt i;
  HighsInt j;

  for (i = 0; i < n / 2; i++) {
    j = a[i];
    a[i] = a[n - 1 - i];
    a[n - 1 - i] = j;
  }

  return;
}
//****************************************************************************80

static void Highs_rcm_level_set(HighsInt root, HighsInt adj_num,
                                const HighsInt adj_row[], const HighsInt adj[],
                                HighsInt mask[], HighsInt* level_num,
                                HighsInt level_row[], HighsInt level[],
                                HighsInt node_num)

//****************************************************************************80
//
//  Purpose:
//
//    level_set() generates the connected level structure rooted at a given
//    node.
//
//  Discussion:
//
//    Only nodes for which MASK is nonzero will be considered.
//
//    The root node chosen by the user is assigned level 1, and masked.
//    All (unmasked) nodes reachable from a node in level 1 are
//    assigned level 2 and masked.  The process continues until there
//    are no unmasked nodes adjacent to any node in the current level.
//    The number of levels may vary between 2 and NODE_NUM.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    This version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Input:
//
//    int ROOT, the node at which the level structure
//    is to be rooted.
//
//    int ADJ_NUM, the number of adjacency entries.
//
//    int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    int MASK[NODE_NUM]: only nodes with nonzero MASK are to be processed.
//
//    int NODE_NUM, the number of nodes.
//
//  Output:
//
//    int MASK[NODE_NUM]: those nodes which were included
//    in the level set have MASK set to 1.
//
//    int *LEVEL_NUM, the number of levels in the level
//    structure.  ROOT is in level 1.  The neighbors of ROOT
//    are in level 2, and so on.
//
//    int LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted
//    level structure.
//
{
  HighsInt i;
  HighsInt iccsze;
  HighsInt j;
  HighsInt jstop;
  HighsInt jstrt;
  HighsInt lbegin;
  HighsInt lvlend;
  HighsInt lvsize;
  HighsInt nbr;
  HighsInt node;

  mask[root - 1] = 0;
  level[0] = root;
  *level_num = 0;
  lvlend = 0;
  iccsze = 1;
  //
  //  LBEGIN is the pointer to the beginning of the current level, and
  //  LVLEND points to the end of this level.
  //
  for (;;) {
    lbegin = lvlend + 1;
    lvlend = iccsze;
    *level_num = *level_num + 1;
    level_row[*level_num - 1] = lbegin;
    //
    //  Generate the next level by finding all the masked neighbors of nodes
    //  in the current level.
    //
    for (i = lbegin; i <= lvlend; i++) {
      node = level[i - 1];
      jstrt = adj_row[node - 1] + offset;
      jstop = adj_row[node] - 1 + offset;

      for (j = jstrt; j <= jstop; j++) {
        nbr = adj[j - 1] + offset;

        if (mask[nbr - 1] != 0) {
          iccsze = iccsze + 1;
          level[iccsze - 1] = nbr;
          mask[nbr - 1] = 0;
        }
      }
    }
    //
    //  Compute the current level width (the number of nodes encountered.)
    //  If it is positive, generate the next level.
    //
    lvsize = iccsze - lvlend;

    if (lvsize <= 0) {
      break;
    }
  }

  level_row[*level_num] = lvlend + 1;
  //
  //  Reset MASK to 1 for the nodes in the level structure.
  //
  for (i = 0; i < iccsze; i++) {
    mask[level[i] - 1] = 1;
  }

  return;
}
//****************************************************************************80

static HighsInt Highs_rcm(HighsInt root, HighsInt adj_num,
                          const HighsInt adj_row[], const HighsInt adj[],
                          HighsInt mask[], HighsInt perm[], HighsInt* iccsze,
                          HighsInt node_num)

//****************************************************************************80
//
//  Purpose:
//
//    rcm() renumbers a connected component by the reverse Cuthill McKee
//    algorithm.
//
//  Discussion:
//
//    The connected component is specified by a node ROOT and a mask.
//    The numbering starts at the root node.
//
//    An outline of the algorithm is as follows:
//
//    X(1) = ROOT.
//
//    for ( I = 1 to N-1)
//      Find all unlabeled neighbors of X(I),
//      assign them the next available labels, in order of increasing degree.
//
//    When done, reverse the ordering.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    09 August 2013
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    This version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Input:
//
//    int ROOT, the node that defines the connected component.
//    It is used as the starting point for the RCM ordering.
//
//    int ADJ_NUM, the number of adjacency entries.
//
//    int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    int MASK[NODE_NUM], a mask for the nodes.  Only
//    those nodes with nonzero input mask values are considered by the
//    routine.
//
//    int NODE_NUM, the number of nodes.
//
//  Output:
//
//    int MASK[NODE_NUM]: The nodes numbered by RCM will have their mask values
//    set to zero.
//
//    int PERM[NODE_NUM], the RCM ordering.
//
//    int *ICCSZE, the size of the connected component
//    that has been numbered.
//
//  Local:
//
//    int DEG[NODE_NUM], a temporary vector used to hold
//    the degree of the nodes in the section graph specified by mask and root.
//
{
  HighsInt* deg;
  HighsInt fnbr;
  HighsInt i;
  HighsInt j;
  HighsInt jstop;
  HighsInt jstrt;
  HighsInt k;
  HighsInt l;
  HighsInt lbegin;
  HighsInt lnbr;
  HighsInt lperm;
  HighsInt lvlend;
  HighsInt nbr;
  HighsInt node;
  //
  //  If node_num out of bounds, something is wrong.
  //
  if (node_num < 1) {
    HIGHS_ORDERING_PRINT("RCM(): Fatal error!\n");
    HIGHS_ORDERING_PRINT("  Unacceptable input value of NODE_NUM\n");
    return 1;
  }
  //
  //  If the root is out of bounds, something is wrong.
  //
  if (root < 1 || node_num < root) {
    HIGHS_ORDERING_PRINT("RCM(): Fatal error!\n");
    HIGHS_ORDERING_PRINT("  Unacceptable input value of ROOT\n");
    return 1;
  }
  //
  //  Allocate memory for the degree array.
  //
  deg = new HighsInt[node_num];
  //
  //  Find the degrees of the nodes in the component specified by MASK and ROOT.
  //
  Highs_rcm_degree(root, adj_num, adj_row, adj, mask, deg, iccsze, perm,
                   node_num);
  //
  //  If the connected component size is less than 1, something is wrong.
  //
  if (*iccsze < 1) {
    HIGHS_ORDERING_PRINT("RCM(): Fatal error!\n");
    HIGHS_ORDERING_PRINT(
        "  Unacceptable connected component size returned from DEGREE\n");
    return 1;
  }
  //
  //  Set the mask value for the root.
  //
  mask[root - 1] = 0;
  //
  //  If the connected component is a singleton, there is no ordering necessary.
  //
  if (*iccsze == 1) {
    delete[] deg;
    return 0;
  }
  //
  //  Carry out the reordering.
  //
  //  LBEGIN and LVLEND point to the beginning and
  //  the end of the current level respectively.
  //
  lvlend = 0;
  lnbr = 1;

  while (lvlend < lnbr) {
    lbegin = lvlend + 1;
    lvlend = lnbr;

    for (i = lbegin; i <= lvlend; i++) {
      //
      //  For each node in the current level...
      //
      node = perm[i - 1];
      jstrt = adj_row[node - 1] + offset;
      jstop = adj_row[node] - 1 + offset;
      //
      //  Find the unnumbered neighbors of NODE.
      //
      //  FNBR and LNBR point to the first and last neighbors
      //  of the current node in PERM.
      //
      fnbr = lnbr + 1;

      for (j = jstrt; j <= jstop; j++) {
        nbr = adj[j - 1] + offset;

        if (mask[nbr - 1] != 0) {
          lnbr = lnbr + 1;
          mask[nbr - 1] = 0;
          perm[lnbr - 1] = nbr;
        }
      }
      //
      //  If no neighbors, skip to next node in this level.
      //
      if (lnbr <= fnbr) {
        continue;
      }
      //
      //  Sort the neighbors of NODE in increasing order by degree.
      //  Linear insertion is used.
      //
      k = fnbr;

      while (k < lnbr) {
        l = k;
        k = k + 1;
        nbr = perm[k - 1];

        while (fnbr < l) {
          lperm = perm[l - 1];

          if (deg[lperm - 1] <= deg[nbr - 1]) {
            break;
          }

          perm[l] = lperm;
          l = l - 1;
        }
        perm[l] = nbr;
      }
    }
  }
  //
  //  We now have the Cuthill-McKee ordering.
  //  Reverse it to get the Reverse Cuthill-McKee ordering.
  //
  Highs_rcm_i4vec_reverse(*iccsze, perm);
  //
  //  Free memory.
  //
  delete[] deg;

  return 0;
}
//****************************************************************************80

static void Highs_rcm_root_find(HighsInt* root, HighsInt adj_num,
                                const HighsInt adj_row[], const HighsInt adj[],
                                HighsInt mask[], HighsInt* level_num,
                                HighsInt level_row[], HighsInt level[],
                                HighsInt node_num)

//****************************************************************************80
//
//  Purpose:
//
//    root_find() finds a pseudo-peripheral node.
//
//  Discussion:
//
//    The diameter of a graph is the maximum distance (number of edges)
//    between any two nodes of the graph.
//
//    The eccentricity of a node is the maximum distance between that
//    node and any other node of the graph.
//
//    A peripheral node is a node whose eccentricity equals the
//    diameter of the graph.
//
//    A pseudo-peripheral node is an approximation to a peripheral node;
//    it may be a peripheral node, but all we know is that we tried our
//    best.
//
//    The routine is given a graph, and seeks pseudo-peripheral nodes,
//    using a modified version of the scheme of Gibbs, Poole and
//    Stockmeyer.  It determines such a node for the section subgraph
//    specified by MASK and ROOT.
//
//    The routine also determines the level structure associated with
//    the given pseudo-peripheral node; that is, how far each node
//    is from the pseudo-peripheral node.  The level structure is
//    returned as a list of nodes LS, and pointers to the beginning
//    of the list of nodes that are at a distance of 0, 1, 2, ...,
//    NODE_NUM-1 from the pseudo-peripheral node.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    This version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//    Norman Gibbs, William Poole, Paul Stockmeyer,
//    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
//    SIAM Journal on Numerical Analysis,
//    Volume 13, pages 236-250, 1976.
//
//    Norman Gibbs,
//    Algorithm 509: A Hybrid Profile Reduction Algorithm,
//    ACM Transactions on Mathematical Software,
//    Volume 2, pages 378-387, 1976.
//
//  Input:
//
//    int *ROOT: a node in the the component of the graph for which a
//    pseudo-peripheral node is sought.
//
//    int ADJ_NUM, the number of adjacency entries.
//
//    int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    int MASK[NODE_NUM], specifies a section subgraph.  Nodes
//    for which MASK is zero are ignored by FNROOT.
//
//    int NODE_NUM, the number of nodes.
//
//  Output:
//
//    int *ROOT: the pseudo-peripheral node obtained.
//
//    int *LEVEL_NUM, is the number of levels in the level structure
//    rooted at the node ROOT.
//
//    int LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
//    level structure array pair containing the level structure found.
//
{
  HighsInt iccsze;
  HighsInt j;
  HighsInt jstrt;
  HighsInt k;
  HighsInt kstop;
  HighsInt kstrt;
  HighsInt level_num2;
  HighsInt mindeg;
  HighsInt nabor;
  HighsInt ndeg;
  HighsInt node;
  //
  //  Determine the level structure rooted at ROOT.
  //
  Highs_rcm_level_set(*root, adj_num, adj_row, adj, mask, level_num, level_row,
                      level, node_num);
  //
  //  Count the number of nodes in this level structure.
  //
  iccsze = level_row[*level_num] - 1;
  //
  //  Extreme case:
  //    A complete graph has a level set of only a single level.
  //    Every node is equally good (or bad).
  //
  if (*level_num == 1) {
    return;
  }
  //
  //  Extreme case:
  //    A "line graph" 0--0--0--0--0 has every node in its only level.
  //    By chance, we've stumbled on the ideal root.
  //
  if (*level_num == iccsze) {
    return;
  }
  //
  //  Pick any node from the last level that has minimum degree
  //  as the starting point to generate a new level set.
  //
  for (;;) {
    mindeg = iccsze;

    jstrt = level_row[*level_num - 1];
    *root = level[jstrt - 1];

    if (jstrt < iccsze) {
      for (j = jstrt; j <= iccsze; j++) {
        node = level[j - 1];
        ndeg = 0;
        kstrt = adj_row[node - 1] + offset;
        kstop = adj_row[node] - 1 + offset;

        for (k = kstrt; k <= kstop; k++) {
          nabor = adj[k - 1] + offset;
          if (0 < mask[nabor - 1]) {
            ndeg = ndeg + 1;
          }
        }

        if (ndeg < mindeg) {
          *root = node;
          mindeg = ndeg;
        }
      }
    }
    //
    //  Generate the rooted level structure associated with this node.
    //
    Highs_rcm_level_set(*root, adj_num, adj_row, adj, mask, &level_num2,
                        level_row, level, node_num);
    //
    //  If the number of levels did not increase, accept the new ROOT.
    //
    if (level_num2 <= *level_num) {
      break;
    }

    *level_num = level_num2;
    //
    //  In the unlikely case that ROOT is one endpoint of a line graph,
    //  we can exit now.
    //
    if (iccsze <= *level_num) {
      break;
    }
  }

  return;
}
//****************************************************************************80

HighsInt Highs_genrcm(HighsInt node_num, HighsInt adj_num,
                      const HighsInt adj_row[], const HighsInt adj[],
                      HighsInt perm[])

//****************************************************************************80
//
//  Purpose:
//
//    genrcm() finds the reverse Cuthill-Mckee ordering for a general graph.
//
//  Discussion:
//
//    For each connected component in the graph, the routine obtains
//    an ordering by calling RCM.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    This version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Input:
//
//    int NODE_NUM, the number of nodes.
//
//    int ADJ_NUM, the number of adjacency entries.
//
//    int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    int  ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    NB: adj_row and adj should be provided using zero-based numbering.
//
//  Output:
//
//    int  PERM[NODE_NUM], the RCM ordering.
//
//    NB: perm is returned using zero-based numbering.
//
//  Local:
//
//    int  LEVEL_ROW[NODE_NUM+1], the index vector for a level
//    structure.  The level structure is stored in the currently unused
//    spaces in the permutation vector PERM.
//
//    int MASK[NODE_NUM], marks variables that have been numbered.
//
{
  HighsInt i;
  HighsInt iccsze;
  HighsInt level_num;
  HighsInt* level_row;
  HighsInt* mask;
  HighsInt num;
  HighsInt root;

  level_row = new HighsInt[node_num + 1];
  mask = new HighsInt[node_num];

  for (i = 0; i < node_num; i++) {
    mask[i] = 1;
  }

  num = 1;

  for (i = 0; i < node_num; i++) {
    //
    //  For each masked connected component...
    //
    if (mask[i] != 0) {
      root = i + 1;
      //
      //  Find a pseudo-peripheral node ROOT.  The level structure found by
      //  ROOT_FIND is stored starting at PERM(NUM).
      //
      Highs_rcm_root_find(&root, adj_num, adj_row, adj, mask, &level_num,
                          level_row, perm + num - 1, node_num);
      //
      //  RCM orders the component using ROOT as the starting node.
      //
      if (Highs_rcm(root, adj_num, adj_row, adj, mask, perm + num - 1, &iccsze,
                    node_num)) {
        delete[] level_row;
        delete[] mask;
        return 1;
      }

      num = num + iccsze;
      //
      //  We can stop once every node is in one of the connected components.
      //
      if (node_num < num) {
        break;
      }
    }
  }

  // convert to zero-based numbering
  for (HighsInt i = 0; i < node_num; ++i) {
    perm[i] -= offset;
  }

  delete[] level_row;
  delete[] mask;

  return 0;
}
//****************************************************************************80
