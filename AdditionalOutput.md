Additional  outputs:
--------------------------

* ```distance_from_leader``` - an array of distances from each cell to the first cell that died in the same experiment.

* ```neigbors_list``` - an array containing a list of neighbors for each cell. The neighbors for cell whose index is i can be found in cell i of the array.

* ```neigbors_list2``` - an array containing a list of first and second degree neighbors for each cell. The neighbors for cell whose index is i can be found in cell i of the array.

* ```neigbors_list3``` - an array containing a list of first, second and third degree neighbors for each cell. The neighbors for cell whose index is i can be found in cell i of the array.

* ```local_leader_flag``` - an array containing the local ‘leader’ for each cell, where a local leader defines the cell that died first out of all of its neighbors.

* ```death_wave_from_local_leader``` – an array that contains the distance of cells from their closest local leader (see definition above).
