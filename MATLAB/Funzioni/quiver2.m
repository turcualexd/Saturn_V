function quiver2(pos, color)

a = annotation('arrow', 'HeadStyle', 'plain', 'HeadLength', 5, 'HeadWidth', 5, 'LineWidth', 1, 'Color', color);
set(a, 'parent', gca);
set(a, 'position', pos);