from Structure2D import Structure2D
# import sympy as sp

import matplotlib.pyplot as plt
import os  # to save plot only

s = Structure2D()

s.add_support(0, 0)

s.add_member_coordinates(s.supports[0].x, s.supports[0].y, s.supports[0].x, 2.1)
s.add_member_length_angle(s.nodes[1].x, s.nodes[1].y, 2.8, 0)
s.add_member_coordinates(s.nodes[0].x, s.nodes[0].y, s.nodes[2].x, s.nodes[2].y)
s.add_member_coordinates(s.nodes[1].x, s.nodes[1].y, 6.4, -2.7)
s.add_member_coordinates(s.nodes[2].x, s.nodes[2].y, 6.4, 2.1)
s.add_member_coordinates(s.nodes[2].x, s.nodes[2].y,2.8+1.575,0)
s.add_member_coordinates(s.nodes[5].x, s.nodes[5].y, 6.4, -2.7)
s.add_member_coordinates(s.nodes[4].x, s.nodes[4].y, s.nodes[5].x, s.nodes[5].y)
s.add_member_coordinates(s.nodes[4].x, s.nodes[4].y, s.nodes[3].x, s.nodes[3].y)
s.add_support(6.4,-2.7,"pin")


for node in s.nodes:
    s.add_hinge(node.x, node.y)

s.add_point_load_global(0,2.1, 1.5, 0,draw_at_head=False)
s.add_point_load_local('n2', None, 1)
s.add_point_load_local('n4', None, 1)
s.add_point_load_local('m3',s.members[3].length/2, 1, s.members[3].angle_deg + 90,draw_at_head=False)


##### plot the structure ####
plt.figure(figsize=(11, 11))
s.plot(
    draw_size=2,
    show_axis=True,
    show_properties=False,
    show_member_labels=True,  # show labels under each member
    draw_all_nodes=True,
    show_node_labels=True,  # show node labels on top of each node
    show_node_labels_types=False,  # show node type
    show_length_bars=True,
    length_bar_exact_values=False,  # True WORKS WITH SYMPY :) , False is decimal values
)

# save the plot
current_directory = os.path.dirname(os.path.abspath(__file__))
save_path = os.path.join(current_directory, "example2.svg")

plt.savefig(save_path, bbox_inches="tight", pad_inches=0.1)

