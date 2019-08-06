function val = generate_contact_matrix(b, num_grps)
% Allocate memories
val = zeros(num_grps);

% Assing values
val(1, :) = [b(1) * ones(1, 19), b(2) * ones(1, 11)];
val(2, 2:end) = [b(1) * ones(1, 18), b(2) * ones(1, 11)];
val(3, 3:end) = [b(2) * ones(1, 17), b(3) * ones(1, 11)];
val(4, 4:end) = [b(2) * ones(1, 16), b(3) * ones(1, 11)];
val(5, 5:end) = b(3) * ones(1, 26);
val(6, 6:end) = b(3) * ones(1, 25);
val(7, 7:end) = b(3) * ones(1, 24);
val(8, 8:end) = [b(4) * ones(1, 7), b(3) * ones(1, 16)];
val(9, 9:end) = [b(4) * ones(1, 6), b(3) * ones(1, 16)];
val(10, 10:end) = [b(5) * ones(1, 5), b(6) * ones(1, 5), b(3) * ones(1, 11)];
val(11, 11:end) = [b(5) * ones(1, 4), b(6) * ones(1, 5), b(3) * ones(1, 11)];
val(12, 12:end) = [b(5) * ones(1, 3), b(6) * ones(1, 5), b(3) * ones(1, 11)];
val(13, 13:end) = [b(5) * ones(1, 2), b(6) * ones(1, 5), b(3) * ones(1, 11)];
val(14, 14:end) = [b(5), b(6) * ones(1, 5), b(3) * ones(1, 11)];
val(15, 15:end) = [b(6) * ones(1, 5), b(3) * ones(1, 11)];
val(16, 16:end) = [b(6) * ones(1, 4), b(7) * ones(1, 4), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(17, 17:end) = [b(6) * ones(1, 3), b(7) * ones(1, 4), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(18, 18:end) = [b(6) * ones(1, 2), b(7) * ones(1, 4), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(19, 19:end) = [b(6), b(7) * ones(1, 4), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(20, 20:end) = [b(7) * ones(1, 4), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(21, 21:end) = [b(7) * ones(1, 3), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(22, 22:end) = [b(7) * ones(1, 2), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(23, 23:end) = [b(7), b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(24, 24:end) = [b(8) * ones(1, 4), b(9) * ones(1, 3)];
val(25, 25:end) = [b(8) * ones(1, 3), b(9) * ones(1, 3)];
val(26, 26:end) = [b(8) * ones(1, 2), b(9) * ones(1, 3)];
val(27, 27:end) = [b(8), b(9) * ones(1, 3)];
val(28, 28:end) = b(9) * ones(1, 3);
val(29, 29:end) = b(9) * ones(1, 2);
val(30, 30:end) = b(9);



