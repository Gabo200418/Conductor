from flask import Flask, render_template, request 
import numpy as np
import math

app = Flask(__name__)#instancia de Flask
@app.route('/') #ruta para la página principal
def index():
    return render_template ('index.html')

@app.route('/Urp_Ucw', methods=['POST']) # Definimos una ruta para manejar el proceso de multiplicación
def Urp_Ucw():
    if request.method == 'POST':
        #Parámetros del sistema
        MVA = float(request.form['num1']) #1
        vN = float(request.form['num2']) #2
        Cn_A = np.round(MVA*1000/(np.sqrt(3)*vN),2) #3
        #Cst = Cn_A #4
        f_d = float(request.form['num23']) #2
        Cd_A = Cn_A + Cn_A * f_d/100
        Cd_A_redondeo = round(Cd_A, 2)
        A_km_snm = float(request.form['num3']) #6
        if A_km_snm < 0:
            resp1 = "No se puede ingresar un número negativo. Por favor, ingrese un valor válido."
        else:
            resp1 = ""

        # Resistencia DC, perdida de calor por convección y radiación
        # Obtener el valor del formulario
        An_c_m = float(request.form['num4'])
        # Verificar si el número es negativo
        if An_c_m < 0:
    
            resp = "No se puede ingresar un número negativo. Por favor, ingrese un valor válido."            
        else:            
            resp = ""

        Al_c_m = float(request.form['num5'])
        if Al_c_m < 0:
            res = "No se puede ingresar un número negativo. Por favor, ingrese un valor válido."
        else:
            res = ""
        IACS = float(request.form['num6'])
        E_ex = float(request.form['num7'])
        Ac =  Al_c_m * An_c_m
        Ta = float(request.form['num8'])
        T2 = float(request.form['num9'])
        R = (1.724 * 10**-6)/ ((IACS * Ac))*(1+((0.00393 * IACS /100)*(T2 - 20)))
        Av = 2 * Al_c_m
        Ah = 2 * An_c_m
        qch_fcv = 3.05 * np.sqrt(1/Al_c_m) * Av * (T2 - Ta)
        qch_fch = 3.05 * np.sqrt(1/An_c_m) * Ah * (T2 - Ta)
        qcv_ncv = 1.38 * Av * ((T2 - Ta) ** 1.25) / (Al_c_m ** 0.25)
        qch_nch_uf = 1.38 * An_c_m * (T2 - Ta) ** 1.25 / An_c_m ** 0.25
        qch_nch_df = 0.69 * An_c_m * (T2 - Ta) ** 1.25 / An_c_m ** 0.25
        qr_fc = qch_fch + qch_fcv
        qc_nc_wm = qch_nch_df + qcv_ncv + qch_nch_uf
        qr_rl_wm = 5.6697 * 10**-8 * E_ex *((Av + Ah))*((T2 + 273)**4-(Ta +273)**4)
        sin_texto = 5.6697*10**-8*0.5*2*(0.0127+0.1524)*(353**4-312**4)


        ## Convertir grados a radianes
        #radianes_B44 = math.radians(Hc)
        #radianes_B45_B46 = math.radians(Zc - Z1)


        #Ganacia de Calor solar
        Lns_gr = float(request.form['num10'])
        Lns_rad = math.radians(Lns_gr)
        Hr_dia = float(request.form['num11'])
        He = A_km_snm * 1000
        w = math.radians((Hr_dia - 12) * 15)
        N = float(request.form['num12'])
        lds_rad = math.radians(23.46 * math.sin(math.radians((284 + N) * 360 / 365)))
        Vas_gr = math.degrees(math.atan(math.sin(w) / (math.sin(Lns_rad) * math.cos(w) - math.cos(Lns_rad) * math.tan(lds_rad))))
        
        if math.degrees(w) < 0:
            Cas_gr = 0 if Vas_gr >= 0 else 180
        else:
            Cas_gr = 180 if Vas_gr >= 0 else 360

        Hc = math.degrees(math.asin(math.cos(Lns_rad)*math.cos(lds_rad)*math.cos(w) + math.sin(Lns_rad)*math.sin(lds_rad)))
        Qs_l =   -42.2391 + 63.8044 * Hc + (-1.922 * Hc**2) + 0.0346921 * Hc**3 + (-0.000361118 * Hc**4) + 0.00000194318 * Hc**5 + (-0.00000000407608 * Hc**6)
        Qs_i =  53.1821 + 14.211 * Hc + (0.66138 * Hc**2) + (-0.031658 * Hc**3) + (0.00054654 * Hc**4) + (-0.0000043446 * Hc**5) + (0.000000013236 * Hc**6)
        AApcul_m2 = Al_c_m * math.sin(math.radians(90 - Hc)) + An_c_m * math.sin(math.radians(Hc))
        k = 1 + 0.0001148 * He + (-0.00000001108 * He**2)
        Zc = Cas_gr + Vas_gr
        Z1 = float(request.form['num13'])
        Aeis_rad = math.acos(math.cos(math.radians(Hc)) * math.cos(math.radians(Zc - Z1)))
        Aeis_gr = math.degrees(Aeis_rad)
        qs_crs_w = E_ex * Qs_l * AApcul_m2 * k * math.sin(Aeis_rad)
        sin_texto1 = math.sin(qs_crs_w)

        #Capacidad de Corriente
        F_cep = float(request.form['num14'])
        I_cc_A =  math.sqrt((qr_fc + qr_rl_wm - qs_crs_w) / (R * F_cep))

        #Gradiente de Tensión
        C_cons = 3.41*10**-4
        Ac_mm2 = (An_c_m * 1000) * (Al_c_m * 1000)
        t_df_s = float(request.form['num15'])
        Ti_tic_c = float(request.form['num16'])
        Tf_tfpc_c = float(request.form['num17'])
        icc_ccc_a = C_cons * Ac_mm2 * 1000000 * math.sqrt((1/t_df_s) * math.log10((Tf_tfpc_c - 20 + (25400/IACS)) / (Ti_tic_c - 20 + (25400/IACS))))

        #Gradiente de Tensión
        rc_rec_cm = Al_c_m*100/2
        Eo_ce = float(request.form['num18'])
        m_fic = float(request.form['num19'])
        c_ce = float(request.form['num20'])
        To_tar_c = float(request.form['num21'])
        Da_dra = ((273+To_tar_c)/(273+Ta))*(1-A_km_snm/10)
        
        Ec_gic = m_fic*Eo_ce*Da_dra*(1+c_ce/(math.sqrt(Da_dra*rc_rec_cm)))

        #Gradiente de tensión el conductor
        H_dcct_cm = float(request.form['num22'])
        D_fd_cm = float(request.form['num23'])

        he_dectf_cm = (H_dcct_cm * D_fd_cm) / math.sqrt(4 * (H_dcct_cm**2) + (D_fd_cm**2))
        d_dc_cm = Al_c_m*100
        V_tns_ft_kv = vN / math.sqrt(3)
        v1_tns_ft_kv = V_tns_ft_kv*1.1
        Ea_gtp = v1_tns_ft_kv / ((d_dc_cm/2) * math.log(4 * he_dectf_cm / d_dc_cm))
        Em_gtm = he_dectf_cm*Ea_gtp/(he_dectf_cm-d_dc_cm/2)
        
        if Em_gtm < Ec_gic:
            Em_Ec = "Cumple"
        else:
            Em_Ec = "No cumple"

            resp = "No se puede ingresar un número negativo en num1. Por favor, ingrese un valor válido."
        return render_template('index.html', r1=Cn_A, r3=Cd_A_redondeo,r0=resp1, r00=resp, r01=res, r4=Ac, r5=R, r6=Av, r7=Ah, r8=qch_fcv, r9=qch_fch, r10=qcv_ncv, r11=qch_nch_uf, r12=qch_nch_df, r13=qr_fc, r14=qc_nc_wm, r15=qr_rl_wm, r16=sin_texto,
                              r17=Lns_rad, r18=He, r19=w, r20=lds_rad, r21=Vas_gr, r22=Cas_gr, r23=Qs_l, r24=Qs_i, r25=AApcul_m2, r26=k, r27=Hc, r28=Zc, r29=Aeis_rad, r30=Aeis_gr, r31=qs_crs_w, r32=sin_texto1, r33=I_cc_A,
                              r34=C_cons, r35=Ac_mm2, r36=icc_ccc_a, r37=rc_rec_cm, r38=Da_dra, r39=Ec_gic, r40=he_dectf_cm, r41=d_dc_cm, r42=V_tns_ft_kv, r43=v1_tns_ft_kv, r44=Ea_gtp, r45=Em_gtm, r46=Em_Ec)
    
if __name__ == '__main__':
    app.run(debug=True)     



    