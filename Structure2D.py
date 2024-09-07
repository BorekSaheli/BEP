"""
This module can be used to solve a 2D structure bending problem with
singularity functions in mechanics.
"""

from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import atan2
from sympy.geometry.polygon import deg
from sympy.simplify.simplify import simplify
from sympy.simplify import nsimplify
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.trigonometric import sin, cos
from sympy.core.numbers import pi

# from sympy.plotting import plot,

# from matplotlib.patches import FancyArrow


import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

import os  # to save plot only


class Member:
    def __init__(self, x1, y1, x2, y2, E, I_flex_rigid, A, member_id):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.E = E
        self.I_flex_rigid = I_flex_rigid
        self.A = A
        self.length = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        self.angle = atan2(y2 - y1, x2 - x1)
        self.angle_deg = deg(self.angle)
        self.member_id = member_id

        # forces
        self.external_loads = []


class ExternalLoad:
    def __init__(
        self,
        load_type,
        x1,
        y1,
        x2=None,
        y2=None,
        applied_to=None,
        local_x=None,
        value=1,
        global_angle=90,
        relative_angle=None,
        load_id=None,
    ):
        self.load_type = load_type
        self.x1 = x1
        self.y1 = y1

        # self.x = x1
        # self.y = y1

        # for distributed loads
        self.x2 = x2
        self.y2 = y2

        self.value = value  # kn
        self.global_angle = global_angle
        self.relative_angle = relative_angle

        self.applied_to = applied_to

        self.local_x = local_x  # Position along the member (used for distributed loads or point loads on members)
        self.load_id = load_id


class Node:
    def __init__(self, x, y, node_type, node_id):
        self.x = x
        self.y = y
        self.node_type = node_type
        self.node_id = node_id

        # forces
        self.external_loads = []


class Support:
    def __init__(self, x, y, support_type, global_angle, support_id):
        self.x = x
        self.y = y
        self.support_type = support_type
        self.global_angle = global_angle
        self.support_id = support_id


# _____________________________________________________________________________________________________________________________________


class Structure2D:
    def __init__(self):
        self.members = []
        self.supports = []
        self.nodes = []

    # find xy position at a member based on local x and member id
    # input example "m0"
    def find_xy_at_target(self, target_name, local_x=None):
        # is not very nice but it works
        if target_name[0] == "m":
            member_id_input = int(target_name.replace("m", ""))
            member = self.members[member_id_input]

            x1, y1, x2, y2 = member.x1, member.y1, member.x2, member.y2
            member_length = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            x = x1 + local_x * (x2 - x1) / member_length
            y = y1 + local_x * (y2 - y1) / member_length

            return x, y

        elif target_name[0] == "n":
            node_id_input = int(target_name.replace("n", ""))
            node = self.nodes[node_id_input]

            return node.x, node.y

    def add_point_load_local(self, target, local_x, value, global_angle=90):
        x, y = self.find_xy_at_target(target, local_x)

        self.add_point_load_global(x, y, value, global_angle=global_angle)

    def add_point_load_global(self, x, y, value, global_angle=90):
        if value == 0:
            raise ValueError("The load value cannot be zero")

        # Check if the load is applied at an existing node
        for node in self.nodes:
            if node.x == x and node.y == y:
                load_id = len(node.external_loads)
                external_load = ExternalLoad(
                    "point",
                    x,
                    y,
                    applied_to=f"n{node.node_id}",
                    local_x=1,
                    value=value,
                    global_angle=global_angle,
                    load_id=load_id,
                )
                node.external_loads.append(external_load)
                # print(f"Load of {value}kN applied at node {node.node_id} at ({float(x)}, {float(y)})")
                return

        # Check if the load is applied on an existing member
        for member in self.members:
            x1, y1, x2, y2 = member.x1, member.y1, member.x2, member.y2

            # Handle vertical members
            if simplify(x2 - x1) == 0 and simplify(x1 - x) == 0:
                # print("vertical member",simplify(x1-x), member.member_id)

                if min(y1, y2) < y < max(y1, y2):
                    load_id = len(member.external_loads)
                    local_x = Abs(y - y1)
                    external_load = ExternalLoad(
                        "point",
                        x,
                        y,
                        applied_to=f"m{member.member_id}",
                        local_x=local_x,
                        value=value,
                        global_angle=global_angle,
                        load_id=load_id,
                    )
                    member.external_loads.append(external_load)
                    # print(f"Load of {value}kN applied on vertical member {member.member_id} at ({float(x)}, {float(y)})")
                    return

            elif simplify(x2 - x1) != 0:
                a = (y2 - y1) / (x2 - x1)
                b = y1 - a * x1
                y_computed = (a * x) + b

                # pfff this is slow af but its finally evaluating correctly
                # print(simplify(nsimplify(line_eq) - nsimplify(y)) , member.member_id)

                if simplify(nsimplify(y_computed) - nsimplify(y)) == 0:
                    # print("on the member", member.member_id,member.x1,member.y1,member.x2,member.y2)
                    # to ensure you are on the line inbetween the start and end point
                    if (
                        x >= min(x1, x2)
                        and x <= max(x1, x2)
                        and y >= min(y1, y2)
                        and y <= max(y1, y2)
                    ):
                        load_id = len(member.external_loads)
                        local_x = sqrt((x - x1) ** 2 + (y - y1) ** 2)
                        external_load = ExternalLoad(
                            "point",
                            x,
                            y,
                            applied_to=f"m{member.member_id}",
                            local_x=local_x,
                            value=value,
                            global_angle=global_angle,
                            load_id=load_id,
                        )
                        member.external_loads.append(external_load)
                        # print(f"Load of {value}kN applied on member {member.member_id} at ({float(x)}, {float(y)})")
                        return

    # _____________________________________________________________________________________________________________________________________
    ### MEMBER FUNCTIONS ###
    def add_member_length_angle(
        self, start_x, start_y, length, angle, E=1, I_flex_rigid=1, A=1
    ):
        angle = angle / 180 * pi
        x1 = start_x
        y1 = start_y
        x2 = x1 + length * cos(angle)
        y2 = y1 + length * sin(angle)
        return self.add_member_coordinates(x1, y1, x2, y2, E, I_flex_rigid, A)

    def add_member_coordinates(self, x1, y1, x2, y2, E=1, I_flex_rigid=1, A=1):
        member_id = len(self.members)
        member = Member(x1, y1, x2, y2, E, I_flex_rigid, A, member_id)
        self.members.append(member)

        # Check start node
        self.add_or_update_node(x1, y1, "fixed", overwrite=False)
        # Check end node
        self.add_or_update_node(x2, y2, "fixed", overwrite=False)
        return member

    # _____________________________________________________________________________________________________________________________________
    ### NODE FUNCTIONS ###
    # Check if the node already exists or update type
    def add_or_update_node(self, x, y, new_node_type, overwrite=True):
        # Update the node type if it already exists

        for node in self.nodes:
            if node.x == x and node.y == y:
                if overwrite:
                    node.node_type = new_node_type
                return node

        # Add a new node if it does not exist
        node_id = len(self.nodes)
        new_node = Node(x, y, new_node_type, node_id)
        self.nodes.append(new_node)
        return new_node

    # _____________________________________________________________________________________________________________________________________
    ### SUPPORT FUNCTIONS ###
    def add_support(self, x, y, support_type="pin", global_angle=0):
        global_angle = np.deg2rad(global_angle)
        support_id = len(self.supports)
        support = Support(
            x,
            y,
            support_type=support_type,
            global_angle=global_angle,
            support_id=support_id,
        )
        self.supports.append(support)

        # Determine the new node type based on the support type
        ### NEED MORE STUFF ADDED HERE
        if support_type == "pin":
            new_node_type = "hinge_support"
        elif support_type == "fixed":
            new_node_type = "fixed_support"
        else:
            pass

        # Add or update the node at the support location
        self.add_or_update_node(x, y, new_node_type, overwrite=True)

        return support

    # _____________________________________________________________________________________________________________________________________
    ### HINGE FUNCTIONS ###
    def add_hinge(self, x, y):
        node_id = len(self.nodes)
        hinge = Node(x, y, "hinge", node_id)

        # Update node
        self.add_or_update_node(x, y, "hinge", overwrite=True)
        return hinge

    # _____________________________________________________________________________________________________________________________________
    ############################################# PLOTTING FUNCTIONS ###########################3#############################
    def plot(
        self,
        show_member_labels=False,
        show_properties=False,
        show_forces=True,
        draw_size=1,
        show_axis=False,
        draw_all_nodes=False,
        show_node_labels=False,
        show_node_labels_types=False,
        show_length_bars=True,
        length_bar_exact_values=True,
    ):
        # _____________________________________________________________________________________________________________________________________
        #### length side bars ###
        def get_length_ticks():
            x_values_for_length_lines = []
            y_values_for_length_lines = []

            for node in self.nodes:
                x_values_for_length_lines.append(node.x)
                y_values_for_length_lines.append(node.y)

            # get the x and y values of external loads
            for member in self.members:
                for load in member.external_loads:
                    x_values_for_length_lines.append(load.x1)
                    y_values_for_length_lines.append(load.y1)

            simplified_x = [nsimplify(expr) for expr in x_values_for_length_lines]
            simplified_y = [nsimplify(expr) for expr in y_values_for_length_lines]

            unique_x = list(dict.fromkeys(simplified_x))
            unique_y = list(dict.fromkeys(simplified_y))

            x_length_ticks = sorted(unique_x)
            y_length_ticks = sorted(unique_y)

            return x_length_ticks, y_length_ticks

        # Draw forces
        def draw_force_vector(ax, x, y, size=1, length=2, angle=90, color="darkred"):
            angle = float(angle - 180)
            dx = np.cos(np.radians(angle)) * length
            dy = np.sin(np.radians(angle)) * length

            arrow_patch = patches.FancyArrow(
                x - dx,
                y - dy,
                dx,
                dy,
                length_includes_head=True,
                width=0.01 * size,
                head_width=0.08 * size,
                head_length=0.16 * size,
                color=color,
            )
            ax.add_patch(arrow_patch).zorder = 100

        # Draw nodes
        def draw_nodes(ax, x, y, node_type, size, edgecolor="black", facecolor="black"):
            if node_type == "hinge" or node_type == "hinge_support":
                size = size * 1 / 50
                circle = patches.Circle(
                    (x, y), size, edgecolor=edgecolor, facecolor="white"
                )
                ax.add_patch(circle).zorder = 11
            elif node_type == "fixed" or node_type == "fixed_support":
                size = size * 1 / 25
                rectangle = patches.Rectangle(
                    (x - size / 2, y - size / 2),
                    size,
                    size,
                    edgecolor=edgecolor,
                    facecolor="white",
                )
                ax.add_patch(rectangle).zorder = 11

        # Draw supports
        def draw_support_pin(
            ax, x, y, global_angle, size, edgecolor="black", facecolor="grey"
        ):
            size = size * 1 / 5

            triangle = np.array(
                [
                    [x, y],
                    [x - np.cos(1) * size, y - np.sin(1) * size],
                    [x + np.cos(1) * size, y - np.sin(1) * size],
                ]
            )

            rotation_matrix = np.array(
                [
                    [np.cos(global_angle), -np.sin(global_angle)],
                    [np.sin(global_angle), np.cos(global_angle)],
                ]
            )

            triangle = np.dot(
                triangle - np.array([x, y]), rotation_matrix.T
            ) + np.array([x, y])

            triangle_patch = patches.Polygon(
                triangle, closed=True, edgecolor=edgecolor, facecolor=facecolor
            )

            circle = patches.Circle(
                (x, y), size * 0.1, edgecolor=edgecolor, facecolor="white"
            )

            ax.add_patch(triangle_patch).zorder = 10
            ax.add_patch(circle).zorder = 10

        # _____________________________________________________________________________________________________________________________________
        #### PLOT ####
        # Plot members
        def get_grid_spacing():
            grid_spacing = plt.gca().get_xticks()
            return grid_spacing[1] - grid_spacing[0]

        def plot_member(member):
            x = [float(member.x1), float(member.x2)]
            y = [float(member.y1), float(member.y2)]
            plt.plot(x, y, "k-")

            mid_x = (x[0] + x[1]) / 2
            mid_y = (y[0] + y[1]) / 2
            angle = float(member.angle)

            return mid_x, mid_y, angle

        def adjust_angle_of_text_rad(angle):
            flip_sign = 1
            if (np.pi / 2 < angle < 3 * np.pi / 2) or (
                -3 * np.pi / 2 < angle < -np.pi / 2
            ):
                angle += np.pi
                flip_sign = -1

            return angle, flip_sign

        def plot_members(show_member_labels, show_properties, show_forces, draw_size):
            for member in self.members:
                plot_member(member)
            # this needs to be done after the members are plotted
            grid_spacing = get_grid_spacing()

            # spacing need to be known before plotting the text
            for member in self.members:
                mid_x, mid_y, angle = plot_member(member)
                member_id = member.member_id

                ## more DRY code
                member_id_text = ""

                if show_member_labels and show_properties:
                    member_id_text = (
                        f"$m_{member_id}$ ({member.E},{member.I_flex_rigid},{member.A})"
                    )
                elif show_member_labels:
                    member_id_text = f"$m_{member_id}$"
                elif show_properties:
                    member_id_text = f"({member.E},{member.I_flex_rigid},{member.A})"

                if member_id_text:
                    angle_adjusted_for_text, flip_sign = adjust_angle_of_text_rad(angle)
                    text_x = mid_x + (flip_sign * grid_spacing * 0.10 * np.sin(angle))
                    text_y = mid_y - (flip_sign * grid_spacing * 0.10 * np.cos(angle))

                    plt.text(
                        text_x,
                        text_y,
                        member_id_text,
                        ha="center",
                        va="center",
                        rotation=np.degrees(angle_adjusted_for_text),
                    )

                if show_forces:
                    for load in member.external_loads:
                        x = load.x1
                        y = load.y1
                        value = load.value
                        angle = load.global_angle
                        draw_force_vector(
                            plt.gca(), x, y, size=draw_size, length=value, angle=angle
                        )

        plot_members(
            show_member_labels=show_member_labels,
            show_properties=show_properties,
            show_forces=show_forces,
            draw_size=draw_size,
        )

        # Plot supports
        def plot_supports(draw_size, show_node_labels, show_node_labels_types):
            for support in self.supports:
                x = float(support.x)
                y = float(support.y)
                support_type = support.support_type
                global_angle = float(support.global_angle)

                ### MORE TYPES NEED TO BE ADDED
                if support_type == "pin":
                    draw_support_pin(
                        plt.gca(), x, y, global_angle=global_angle, size=draw_size
                    )
                else:
                    pass

        plot_supports(
            draw_size=draw_size,
            show_node_labels=show_node_labels,
            show_node_labels_types=show_node_labels_types,
        )

        grid_spacing = plt.gca().get_xticks()
        grid_spacing = grid_spacing[1] - grid_spacing[0]
        for node in self.nodes:
            x = float(node.x)
            y = float(node.y)
            node_type = node.node_type

            if node_type == "hinge" or node_type == "support_hinge":
                draw_nodes(
                    plt.gca(), x, y, node_type, size=draw_size, edgecolor="black"
                )

            if draw_all_nodes:
                draw_nodes(
                    plt.gca(), x, y, node_type, size=draw_size, edgecolor="black"
                )

            ## more DRY code
            node_id_text = ""

            if show_node_labels and show_node_labels_types:
                node_id_text = f"$n_{node.node_id}$ $({node.node_type})$"

            elif show_node_labels and not show_node_labels_types:
                node_id_text = f"$n_{node.node_id}$"

            elif show_node_labels_types and not show_node_labels:
                node_id_text = f"$({node.node_type})$"

            if node_id_text:
                plt.text(
                    x,
                    y + grid_spacing / 8,
                    node_id_text,
                    ha="center",
                    va="bottom",
                    color="blue",
                    zorder=12,
                    bbox=dict(facecolor="lightgrey", alpha=0.50, edgecolor="none"),
                )

            if show_forces:
                for load in node.external_loads:
                    x = load.x1
                    y = load.y1
                    value = load.value
                    angle = load.global_angle
                    draw_force_vector(
                        plt.gca(), x, y, size=draw_size, length=value, angle=angle
                    )

        # Draw forces
        if show_length_bars:
            x_length_ticks, y_length_ticks = get_length_ticks()
            x_min = round(float(y_length_ticks[0]), 0)
            y_min = round(float(x_length_ticks[-1]), 0)

            # length lines plot
            # along x
            plt.plot(
                x_length_ticks,
                (x_min - grid_spacing) * np.ones(len(x_length_ticks)),
                "k",
                marker="|",
            )
            # along y
            plt.plot(
                (y_min + grid_spacing) * np.ones(len(y_length_ticks)),
                y_length_ticks,
                "k",
                marker="_",
            )

            # lenth per segment text
            for i in range(0, len(x_length_ticks) - 1):
                loc = (x_length_ticks[i] + x_length_ticks[i + 1]) / 2
                segment_lenth = Abs(simplify(x_length_ticks[i + 1] - x_length_ticks[i]))

                if length_bar_exact_values:
                    # print(float(loc) - y_length_ticks[0])
                    plt.text(
                        float(loc),
                        x_min - grid_spacing * 1.25,
                        f"${segment_lenth}$",
                        ha="center",
                        va="bottom",
                        color="black",
                        fontsize=8,
                    )
                else:
                    plt.text(
                        float(loc),
                        x_min - grid_spacing * 1.25,
                        f"${round(float(segment_lenth),2)}$",
                        ha="center",
                        va="bottom",
                        color="black",
                        fontsize=8,
                    )

            for i in range(0, len(y_length_ticks) - 1):
                loc = (y_length_ticks[i] + y_length_ticks[i + 1]) / 2
                segment_lenth = Abs(simplify(y_length_ticks[i + 1] - y_length_ticks[i]))

                if length_bar_exact_values:
                    # get  grid ticks

                    plt.text(
                        y_min + grid_spacing * 1.25,
                        float(loc),
                        f"${segment_lenth}$",
                        ha="center",
                        va="center",
                        color="black",
                        fontsize=8,
                    )
                else:
                    plt.text(
                        y_min + grid_spacing * 1.25,  # y_min + 1.5,
                        float(loc),
                        f"${round(float(segment_lenth),2)}$",
                        ha="center",
                        va="bottom",
                        color="black",
                        fontsize=8,
                    )
        # _____________________________________________________________________________________________________________________________________
        #### PLOT STYLE ####
        if not show_axis:
            plt.tick_params(
                axis="both",
                which="both",
                bottom=False,
                top=False,
                labelbottom=False,
                right=False,
                left=False,
                labelleft=False,
            )

        ax = plt.gca()
        ax.set_aspect("equal", adjustable="datalim")

        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        plt.grid()


# _____________________________________________________________________________________________________________________________________
#######################################################################################################################################
# _____________________________________________________________________________________________________________________________________
########### EXAMPLE ###################################################################################################################

s = Structure2D()


s.add_support(0, 0)

s.add_member_coordinates(s.supports[0].x, s.supports[0].y, -10, 15)
s.add_member_coordinates(s.members[0].x2, s.members[0].y2, 5, 3.75)
s.add_member_coordinates(
    s.supports[0].x, s.supports[0].y, s.members[1].x2, s.members[1].y2
)

s.add_support(
    10,
    0,
    "pin",
)
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

s.add_point_load_global(20, 15, 4, 180)
s.add_point_load_local("m1", s.members[1].length / 2, 3, s.members[1].angle_deg + 90)
s.add_point_load_local("n2", None, 3, -90)
s.add_point_load_local("m2", 1, 3, 0)

y = 5
a = (s.nodes[6].y - s.nodes[5].y) / (s.nodes[6].x - s.nodes[5].x)
b = s.nodes[6].y - (a * s.nodes[6].x)
x = (y - b) / a

s.add_point_load_global(x, 5, 3, 90)
s.add_point_load_global(x, 5, 3, 90)


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
save_path = os.path.join(current_directory, "structure.svg")


plt.savefig(save_path, bbox_inches="tight", pad_inches=0.1)
