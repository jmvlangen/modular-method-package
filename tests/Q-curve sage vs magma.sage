### 2-isogeny parameter 2 - level 512 = 2^9
E = Qcurve_with_2_isogeny(2)
E = E.decomposable_twist()
E = E.twist(QQ(3))
# Even the playing field:
E.decomposition_field();
E.splitting_character('conjugacy');
E.twist_character('conjugacy');
E.splitting_image_field('conjugacy');
%time f_m, twists = E.newform(algorithm='magma') # 1.17 s
%time f, twists = E.newform(algorithm='sage') # 0.357 s

### 2-isogeny parameter 3 - level 1536 = 2^9 * 3
E = Qcurve_with_2_isogeny(3)
E = E.decomposable_twist()
%time f_m, twists = E.newform(algorithm='magma') # 7.84 s
%time f, twists = E.newform(algorithm='sage') # 1449 s (24 min)

### 2-isogeny parameter -1 - level 512 = 2^9
E = Qcurve_with_2_isogeny(-1)
E = E.decomposable_twist()
E = E.twist(QQ(3))
# Even the playing field:
E.decomposition_field();
E.splitting_character('conjugacy');
E.twist_character('conjugacy');
E.splitting_image_field('conjugacy');
%time f_m, twists = E.newform(algorithm='magma') # 1.43 s
%time f, twists = E.newform(algorithm='sage') # 0.382 s

### 2-isogeny parameter 5 - level 800 = 2^5 * 5^2
E = Qcurve_with_2_isogeny(5)
E = E.decomposable_twist()
E = E.twist(QQ(3))
# Even the playing field:
E.decomposition_field();
E.splitting_character('conjugacy');
E.twist_character('conjugacy');
E.splitting_image_field('conjugacy');
%time f_m, twists = E.newform(algorithm='magma') # 9.65 s
%time f, twists = E.newform(algorithm='sage') # 11.1 s

### 2-isogeny parameter -4 - level 1280 = 2^8 * 5
E = Qcurve_with_2_isogeny(-4)
E = E.decomposable_twist()
E = E.twist(QQ(3))
# Even the playing field:
E.decomposition_field();
E.splitting_character('conjugacy');
E.twist_character('conjugacy');
E.splitting_image_field('conjugacy');
%time f_m, twists = E.newform(algorithm='magma') # 2.42 s
%time f, twists = E.newform(algorithm='sage') # 14.5 s
