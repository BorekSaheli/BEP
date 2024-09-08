from Structure2D import Structure2D
# import sympy as sp

import matplotlib.pyplot as plt
import os  # to save plot only

s = Structure2D()

s.add_support(0, 0)
s.add_member_coordinates(s.supports[0].x, s.supports[0].y, -10, 15)
s.add_member_coordinates(s.members[0].x2, s.members[0].y2, 5, 3.75)
s.add_member_coordinates(s.supports[0].x, s.supports[0].y, s.members[1].x2, s.members[1].y2)
s.add_support(10,0,"pin")
s.add_member_coordinates(s.supports[1].x, s.supports[1].y, s.nodes[2].x, s.nodes[2].y)
s.add_member_coordinates(s.supports[1].x, s.supports[1].y, s.supports[1].x, 7.5)
s.add_member_coordinates(s.nodes[4].x, s.nodes[4].y, s.nodes[2].x, s.nodes[2].y)
s.add_member_coordinates(s.supports[1].x, s.supports[1].y, 12.5, 3.75)
s.add_member_coordinates(s.nodes[5].x, s.nodes[5].y, s.nodes[4].x, s.nodes[4].y)
s.add_member_coordinates(s.nodes[5].x, s.nodes[5].y, 20, 15)
s.add_member_coordinates(s.nodes[6].x, s.nodes[6].y, s.nodes[4].x, s.nodes[4].y)
s.add_support(-10, 15, "pin", 90 + 180)

s.add_hinge(s.nodes[1].x, s.nodes[1].y)
s.add_hinge(s.nodes[2].x, s.nodes[2].y)
s.add_hinge(s.nodes[4].x, s.nodes[4].y)
s.add_hinge(s.nodes[5].x, s.nodes[5].y)
s.add_hinge(s.nodes[6].x, s.nodes[6].y)

s.add_point_load_global(20, 15, 4, 180,draw_at_head=False)
s.add_point_load_local("m1", s.members[1].length / 2, 3, s.members[1].angle_deg + 90)
s.add_point_load_local("n2", None, 3, -90)
s.add_point_load_local("m2", 1, 3, 0)

y = 5
a = (s.nodes[6].y - s.nodes[5].y) / (s.nodes[6].x - s.nodes[5].x)
b = s.nodes[6].y - (a * s.nodes[6].x)
x = (y - b) / a

s.add_point_load_global(x, 5, 3, 90,draw_at_head=False)

##### plot the structure ####
plt.figure(figsize=(11, 11))
s.plot(
    draw_size=8,
    show_axis=True,
    show_properties=False,
    show_member_labels=True,  # show labels under each member
    draw_all_nodes=True,
    show_node_labels=True,  # show node labels on top of each node
    show_node_labels_types=False,  # show node type
    show_length_bars=True,
    length_bar_exact_values=True,  # True WORKS WITH SYMPY :) , False is decimal values
)

# save the plot
current_directory = os.path.dirname(os.path.abspath(__file__))
save_path = os.path.join(current_directory, "example1.svg")

plt.savefig(save_path, bbox_inches="tight", pad_inches=0.1)

print('Memeber 1 external load')
print(f'fx={s.members[8].external_loads[0].global_f_horizontal}')
print(f'fy={s.members[8].external_loads[0].global_f_vertical}')
print(f'angle={s.members[8].external_loads[0].global_angle_rad}')
print(f'angle={s.members[8].external_loads[0].global_angle}')
