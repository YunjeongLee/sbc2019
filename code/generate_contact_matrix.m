function val = generate_contact_matrix(b, num_grps)
% Allocate memories
val = zeros(num_grps);

% Assing values
val(1, :) = [b1 * ones(1, 19), b2 * ones(1, 11)];
val(2, 2:end) = [b1 * one(1, 18), b2 * ones(1, 11)];
val(3, 3:end) = [b2 * ones(1, 17), b3 * ones(1, 11)];
val(4, 4:end) = [b2 * ones(1, 16), b3 * ones(1, 11)];
val(5, 5:end) = b3 * ones(1, 26);
val(6, 6:end) = b3 * ones(1, 25);
val(7, 7:end) = b3 * ones(1, 24);
val(8, 8:end) = [b4 * ones(1, 7), b3 * ones(1, 16)];
val(9, 9:end) = [b4 * ones(1, 6), b3 * ones(1, 16)];
val(10, 10:end) = [b5 * ones(1, 5), b6 * ones(1, 5), b3 * ones(1, 11)];
val(11, 11:end) = [b5 * ones(1, 4), b6 * ones(1, 5), b3 * ones(1, 11)];
val(12, 12:end) = [b5 * ones(1, 3), b6 * ones(1, 5), b3 * ones(1, 11)];
val(13, 13:end) = [b5 * ones(1, 2), b6 * ones(1, 5), b3 * ones(1, 11)];
val(14, 14:end) = [b5, b6 * ones(1, 5), b3 * ones(1, 11)];
val(15, 15:end) = [b6 * ones(1, 5), b3 * ones(1, 11)];
val(16, 16:end) = [b6 * ones(1, 4), b7 * ones(1, 4), b8 * ones(1, 4), b9 * ones(1, 3)];
val(17, 17:end) = [b6 * ones(1, 3), b7 * ones(1, 4), b8 * ones(1, 4), b9 * ones(1, 3)];
val(18, 18:end) = [b6 * ones(1, 2), b7 * ones(1, 4), b8 * ones(1, 4), b9 * ones(1, 3)];
val(19, 19:end) = [b6, b7 * ones(1, 4), b8 * ones(1, 4), b9 * ones(1, 3)];
val(20, 20:end) = [b7 * ones(1, 4), b8 * ones(1, 4), b9 * ones(1, 3)];
val(21, 21:end) = [b7 * ones(1, 3), b8 * ones(1, 4), b9 * ones(1, 3)];
val(22, 22:end) = [b7 * ones(1, 2), b8 * ones(1, 4), b9 * ones(1, 3)];
val(23, 23:end) = [b7, b8 * ones(1, 4), b9 * ones(1, 3)];
val(24, 24:end) = [b8 * ones(1, 4), b9 * ones(1, 3)];
val(25, 25:end) = [b8 * ones(1, 3), b9 * ones(1, 3)];
val(26, 26:end) = [b8 * ones(1, 2), b9 * ones(1, 3)];
val(27, 26:end) = [b8, b9 * ones(1, 3)];
val(28, 26:end) = b9 * ones(1, 3);
val(29, 26:end) = b9 * ones(1, 2);
val(30, 30:end) = b9;



