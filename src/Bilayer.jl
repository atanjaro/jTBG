function generate_extended_unit_cell(m,n,a₁,a₂)
    # total number of points in the super cell
    N_unit_cell_points = n^2 + n*m + m^2

    A1Points = zeros(Int, N_unit_cell_points, 3)
    BondList = zeros(Int, 3 * N_unit_cell_points, 3)

    count = 0
    countB = 0
    for i in -2*(n + m):2*(n + m)
        for j in -2*(n + m):2*(n + m)
            if in_parallelogram([0, 0], [G₁*0.9999999, G₂*0.9999999], TrA(i, j))
                count += 1
                A2Points[count, :] = [i, j, TrA(i, j)]
                countB += 1
                TBondList[countB, :] = [TrA(i, j), TrA(i, j) + Tδ₁]
                countB += 1
                TBondList[countB, :] = [TrA(i, j), TrA(i, j) + Tδ₂]
                countB += 1
                TBondList[countB, :] = [TrA(i, j), TrA(i, j) + Tδ₃]
            end
        end
    end


end