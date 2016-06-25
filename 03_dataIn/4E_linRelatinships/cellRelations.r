# # # # CELL RELATIONSHIP INFORMATION # # # # 


# make a database of lineal relationships for each cell:

linRelations <- list()
cellTypes <- factor(c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", 
                      "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                      "ABalx", "ABarx", "ABplx", "ABprx",
                      "MSx", "MSa", "MSp", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Ex1", "Ex2",
                      "Cx", "Ca", "Cp", "Cx1", "Cx2", "D", "P4"), 
                    levels=c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", 
                             "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                             "ABalx", "ABarx", "ABplx", "ABprx",
                             "MSx", "MSa", "MSp", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Ex1", "Ex2",
                             "Cx", "Ca", "Cp", "Cx1", "Cx2", "D", "P4"))
for(cell in cellTypes){
    linRelations[[cell]] <- list()
}

# DAUGHTERS
linRelations[["P0"]][["daughters"]] <- c("AB", "P1")

linRelations[["P1"]][["daughters"]] <- c("EMS", "P2")
linRelations[["AB"]][["daughters"]] <- c("ABa", "ABp")

linRelations[["ABa"]][["daughters"]] <- c("ABal", "ABar")
linRelations[["ABp"]][["daughters"]] <- c("ABpl", "ABpr")
linRelations[["EMS"]][["daughters"]] <- c("E", "MS")
linRelations[["P2"]][["daughters"]] <- c("C", "P3")

linRelations[["ABal"]][["daughters"]] <- c("ABalx")
linRelations[["ABar"]][["daughters"]] <- c("ABarx")
linRelations[["ABpl"]][["daughters"]] <- c("ABplx")
linRelations[["ABpr"]][["daughters"]] <- c("ABprx")
linRelations[["MS"]][["daughters"]] <- c("MSx1", "MSx2", "MSx")
linRelations[["E"]][["daughters"]] <- c("Ep", "Ea", "Ex")
linRelations[["C"]][["daughters"]] <- c("Cx1", "Cx2", "Cx")
linRelations[["P3"]][["daughters"]] <- c("D", "P4")

# PARENTS
linRelations[["P1"]][["parent"]] <- c("P0")
linRelations[["AB"]][["parent"]] <- c("P0")

linRelations[["ABa"]][["parent"]] <- c("AB")
linRelations[["ABp"]][["parent"]] <- c("AB")
linRelations[["EMS"]][["parent"]] <- c("P1")
linRelations[["P2"]][["parent"]] <- c("P1")

linRelations[["ABal"]][["parent"]] <- c("ABa")
linRelations[["ABar"]][["parent"]] <- c("ABa")
linRelations[["ABpl"]][["parent"]] <- c("ABp")
linRelations[["ABpr"]][["parent"]] <- c("ABp")
linRelations[["MS"]][["parent"]] <- c("EMS")
linRelations[["E"]][["parent"]] <- c("EMS")
linRelations[["C"]][["parent"]] <- c("P2")
linRelations[["P3"]][["parent"]] <- c("P2")

linRelations[["ABalx"]][["parent"]] <- c("ABal")
linRelations[["ABarx"]][["parent"]] <- c("ABar")
linRelations[["ABplx"]][["parent"]] <- c("ABpl")
linRelations[["ABprx"]][["parent"]] <- c("ABpr")
linRelations[["MSx1"]][["parent"]] <- c("MS")
linRelations[["MSx2"]][["parent"]] <- c("MS")
linRelations[["MSx"]][["parent"]] <- c("MS")
linRelations[["Ep"]][["parent"]] <- c("E")
linRelations[["Ea"]][["parent"]] <- c("E")
linRelations[["Ex"]][["parent"]] <- c("E")
linRelations[["Cx1"]][["parent"]] <- c("C")
linRelations[["Cx2"]][["parent"]] <- c("C")
linRelations[["Cx"]][["parent"]] <- c("C")
linRelations[["D"]][["parent"]] <- c("P3")
linRelations[["P4"]][["parent"]] <- c("P3")

# GRAND DAUGHTERS
linRelations[["P0"]][["granddaughters"]] <- c("AB", "P1", "ABa", "ABp", "EMS", "P2")

linRelations[["P1"]][["granddaughters"]] <- c("E", "MS", "C", "P3")
linRelations[["AB"]][["granddaughters"]] <- c("ABal", "ABar", "ABpl", "ABpr")

linRelations[["ABa"]][["granddaughters"]] <- c("ABalx", "ABarx")
linRelations[["ABp"]][["granddaughters"]] <- c("ABplx", "ABprx")
linRelations[["EMS"]][["granddaughters"]] <- c("Ep", "Ea", "Ex", "MSx1", "MSx2", "MSx")
linRelations[["P2"]][["granddaughters"]] <- c("Cx1", "Cx2", "Cx", "D", "P4")

# GRAND PARENTS
linRelations[["ABa"]][["grandparent"]] <- c("P0")
linRelations[["ABp"]][["grandparent"]] <- c("P0")
linRelations[["EMS"]][["grandparent"]] <- c("P0")
linRelations[["P2"]][["grandparent"]] <- c("P0")

linRelations[["ABal"]][["grandparent"]] <- c("AB")
linRelations[["ABar"]][["grandparent"]] <- c("AB")
linRelations[["ABpl"]][["grandparent"]] <- c("AB")
linRelations[["ABpr"]][["grandparent"]] <- c("AB")
linRelations[["MS"]][["grandparent"]] <- c("P1")
linRelations[["E"]][["grandparent"]] <- c("P1")
linRelations[["C"]][["grandparent"]] <- c("P1")
linRelations[["P3"]][["grandparent"]] <- c("P1")

linRelations[["ABalx"]][["grandparent"]] <- c("ABa")
linRelations[["ABarx"]][["grandparent"]] <- c("ABa")
linRelations[["ABplx"]][["grandparent"]] <- c("ABp")
linRelations[["ABprx"]][["grandparent"]] <- c("ABp")
linRelations[["MSx1"]][["grandparent"]] <- c("EMS")
linRelations[["MSx2"]][["grandparent"]] <- c("EMS")
linRelations[["MSx"]][["grandparent"]] <- c("EMS")
linRelations[["Ep"]][["grandparent"]] <- c("EMS")
linRelations[["Ea"]][["grandparent"]] <- c("EMS")
linRelations[["Ex"]][["grandparent"]] <- c("EMS")
linRelations[["Cx1"]][["grandparent"]] <- c("P2")
linRelations[["Cx2"]][["grandparent"]] <- c("P2")
linRelations[["Cx"]][["grandparent"]] <- c("P2")
linRelations[["D"]][["grandparent"]] <- c("P2")
linRelations[["P4"]][["grandparent"]] <- c("P2")

# GREAT GRAND DAUGHTERS
linRelations[["P0"]][["greatgranddaughters"]] <- c("ABal", "ABar", "ABpl", "ABpr", "E", "MS", "C", "P3")

linRelations[["P1"]][["greatgranddaughters"]] <- c("Ep", "Ea", "Ex", "MSx1", "MSx2", "MSx", "Cx1", "Cx2", "Cx", "D", "P4")
linRelations[["AB"]][["greatgranddaughters"]] <- c("ABalx", "ABarx", "ABplx", "ABprx")

# GREAT GRAND PARENTS
linRelations[["ABal"]][["greatgrandparent"]] <- c("P0")
linRelations[["ABar"]][["greatgrandparent"]] <- c("P0")
linRelations[["ABpl"]][["greatgrandparent"]] <- c("P0")
linRelations[["ABpr"]][["greatgrandparent"]] <- c("P0")
linRelations[["MS"]][["greatgrandparent"]] <- c("P0")
linRelations[["E"]][["greatgrandparent"]] <- c("P0")
linRelations[["C"]][["greatgrandparent"]] <- c("P0")
linRelations[["P3"]][["greatgrandparent"]] <- c("P0")

linRelations[["ABalx"]][["greatgrandparent"]] <- c("AB")
linRelations[["ABarx"]][["greatgrandparent"]] <- c("AB")
linRelations[["ABplx"]][["greatgrandparent"]] <- c("AB")
linRelations[["ABprx"]][["greatgrandparent"]] <- c("AB")
linRelations[["MSx1"]][["greatgrandparent"]] <- c("P1")
linRelations[["MSx2"]][["greatgrandparent"]] <- c("P1")
linRelations[["MSx"]][["greatgrandparent"]] <- c("P1")
linRelations[["Ep"]][["greatgrandparent"]] <- c("P1")
linRelations[["Ea"]][["greatgrandparent"]] <- c("P1")
linRelations[["Ex"]][["greatgrandparent"]] <- c("P1")
linRelations[["Cx1"]][["greatgrandparent"]] <- c("P1")
linRelations[["Cx2"]][["greatgrandparent"]] <- c("P1")
linRelations[["Cx"]][["greatgrandparent"]] <- c("P1")
linRelations[["D"]][["greatgrandparent"]] <- c("P1")
linRelations[["P4"]][["greatgrandparent"]] <- c("P1")

# GREAT GREAT GRAND DAUGHTERS
linRelations[["P0"]][["greatgreatgranddaughters"]] <- c("ABalx", "ABarx", "ABplx", "ABprx", "Ep", "Ea", "Ex", "MSx1", "MSx2", "MSx", "Cx1", "Cx2", "Cx", "D", "P4")

# GREAT GREAT GRAND PARENTS
linRelations[["ABalx"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["ABarx"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["ABplx"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["ABprx"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["MSx1"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["MSx2"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["MSx"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["Ep"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["Ea"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["Ex"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["Cx1"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["Cx2"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["Cx"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["D"]][["greatgreatgrandparent"]] <- c("P0")
linRelations[["P4"]][["greatgreatgrandparent"]] <- c("P0")


# SISTERS
linRelations[["P1"]][["sister"]] <- c("AB")
linRelations[["AB"]][["sister"]] <- c("P1")

linRelations[["ABa"]][["sister"]] <- c("ABp")
linRelations[["ABp"]][["sister"]] <- c("ABa")
linRelations[["EMS"]][["sister"]] <- c("P2")
linRelations[["P2"]][["sister"]] <- c("EMS")

linRelations[["ABal"]][["sister"]] <- c("ABar")
linRelations[["ABar"]][["sister"]] <- c("ABal")
linRelations[["ABpl"]][["sister"]] <- c("ABpr")
linRelations[["ABpr"]][["sister"]] <- c("ABpl")
linRelations[["MS"]][["sister"]] <- c("E")
linRelations[["E"]][["sister"]] <- c("MS")
linRelations[["C"]][["sister"]] <- c("P3")
linRelations[["P3"]][["sister"]] <- c("C")

linRelations[["ABalx"]][["sister"]] <- c("ABarx")
linRelations[["ABarx"]][["sister"]] <- c("ABalx")
linRelations[["ABplx"]][["sister"]] <- c("ABprx")
linRelations[["ABprx"]][["sister"]] <- c("ABplx")
linRelations[["MSx1"]][["sister"]] <- c("MSx2")
linRelations[["MSx2"]][["sister"]] <- c("MSx1")
linRelations[["MSx"]][["sister"]] <- c("MSx")
linRelations[["Ep"]][["sister"]] <- c("Ea")
linRelations[["Ea"]][["sister"]] <- c("Ep")
linRelations[["Ex"]][["sister"]] <- c("Ex")
linRelations[["Cx1"]][["sister"]] <- c("Cx2")
linRelations[["Cx2"]][["sister"]] <- c("Cx1")
linRelations[["Cx"]][["sister"]] <- c("Cx")
linRelations[["D"]][["sister"]] <- c("P4")
linRelations[["P4"]][["sister"]] <- c("D")

# COUSINS #

linRelations[["ABa"]][["cousins"]] <- c("EMS", "P2")
linRelations[["ABp"]][["cousins"]] <- c("EMS", "P2")
linRelations[["EMS"]][["cousins"]] <- c("ABa", "ABp")
linRelations[["P2"]][["cousins"]] <- c("ABa", "ABp")

linRelations[["ABal"]][["cousins"]] <- c("ABpl", "ABpr")
linRelations[["ABar"]][["cousins"]] <- c("ABpl", "ABpr")
linRelations[["ABpl"]][["cousins"]] <- c("ABal", "ABar")
linRelations[["ABpr"]][["cousins"]] <- c("ABal", "ABar")
linRelations[["MS"]][["cousins"]] <- c("C", "P3")
linRelations[["E"]][["cousins"]] <- c("C", "P3")
linRelations[["C"]][["cousins"]] <- c("E", "MS")
linRelations[["P3"]][["cousins"]] <- c("E", "MS")

linRelations[["ABalx"]][["cousins"]] <- c("ABplx", "ABprx")
linRelations[["ABarx"]][["cousins"]] <- c("ABplx", "ABprx")
linRelations[["ABplx"]][["cousins"]] <- c("ABalx", "ABarx")
linRelations[["ABprx"]][["cousins"]] <- c("ABalx", "ABarx")
linRelations[["MSx1"]][["cousins"]] <- c("Ep", "Ea")
linRelations[["MSx2"]][["cousins"]] <- c("Ep", "Ea")
linRelations[["MSx"]][["cousins"]] <- c("Ex")
linRelations[["Ep"]][["cousins"]] <- c("MSx1", "MSx2")
linRelations[["Ea"]][["cousins"]] <- c("MSx1", "MSx2")
linRelations[["Ex"]][["cousins"]] <- c("MSx")
linRelations[["Cx1"]][["cousins"]] <- c("D", "P4")
linRelations[["Cx2"]][["cousins"]] <- c("D", "P4")
linRelations[["Cx"]][["cousins"]] <- c("D", "P4")
linRelations[["D"]][["cousins"]] <- c("Cx1", "Cx2", "Cx")
linRelations[["P4"]][["cousins"]] <- c("Cx1", "Cx2", "Cx")


rm(cellTypes, cell)