/*
void test_variant()
{
	for (int i = 0; i < M - 1; i++)
	{
		X[i] = h * i;
	}
	vector <double > X_vec;
	double P01 = 4.0 / 3.0, ro01 = 4.0, e01 = 0.5, U01 = 1.0; // x [0,1]
	double P02 = 0.0002 / 3.0, ro02 = 1.0, e02 = 0.0001, U02 = 0.0; // x (1,2]
	double sigma = 0.6, gamma = 5.0 / 3.0, T = 0.6;
	double time = 0.0, end_time = 0.6, deltaX = 0.0;
	//X = 0.0; //tau - Ð½Ð°Ð´Ð¾ Ð½Ð°Ð¹Ñ‚Ð¸ Ð¿Ñ€Ð°Ð²Ð»ÑŒÐ½Ð¾Ð?
	//Í.Ó.
	for (int i = 0; i < M; i++)
	{
		Q[i] = 0.0;
		if (i <= M / 2)
		{
			U[i] = U01;
			RO[i] = ro01;
			NU[i] = 1 / ro01;
			P[i] = P01;
			C[i] = pow(gamma * P[i] / RO[i], 0.5);
			G[i] = P[i] + Q[i];
			E[i] = e01;

		}
		else
		{
			U[i] = U02;
			RO[i] = ro02;
			NU[i] = 1 / ro02;
			P[i] = P02;
			C[i] = pow(gamma * P[i] / RO[i], 0.5);
			G[i] = P[i] + Q[i];
			E[i] = e02;

		}
	}
	for (int i = 0; i < M - 1; i++)
	{
		deltaS[i] = h * RO[i];
	}


	for (int i = 0; i < M; i++)
	{
		//if (i <= ro01 / deltaS)
		if (i * h <= 1)
		{
			E[i] = e01;
			P[i] = P01;
			P_new[i] = P01;
			C[i] = pow(gamma * P[i] / RO[i], 0.5);
			G[i] = P[i] + Q[i];
		}
		else
		{
			E[i] = e02;
			P[i] = P02;
			C[i] = pow(gamma * P[i] / RO[i], 0.5);
			G[i] = P[i] + Q[i];
		}
	}

	while (time <= T) {
		double min_k = 10;
		double min_uv = 10;
		float C = 0;
		//	//-------------ðàñ÷åò òàó------------
			//for (int i = 0; i < M - 1; i++)
			//{
			//	C = pow(gamma * P[i] / RO[i], 0.5);
			//	tau_k[i] = 10;// 0.1*(deltaS[i] * (sqrt(3 * gamma + 1)*gamma) / (C*(gamma + 1))); //deltaS[i] / C * 4.9;
			//	if (fabs(U[i + 1] - U[i]) >= 0.000001)
			//		tau_uv[i] = 1.0 / (2 * 4 * (fabs(U[i + 1] - U[i])));
			//	else tau_uv[i] = 10;

			//
			//}
		/*	for (int p = 0; p < M - 1; p++) {

				if (tau_k[p] <= min_k) { min_k = tau_k[p]; }
				if (tau_uv[p] <= min_uv) { min_uv = tau_uv[p]; }
				if (min_k <= min_uv) {
					min_k = min_uv;
				}
			}
			if (min_k <= 1.2*tau) tau = min_k;
			//else tau = 1.2*tau;
			//---------------------------------------


		U[0] = 1; U_new[0] = 1;
		U[M - 1] = 0; U_new[M - 1] = 0;

		for (int k = 0; k < M - 1; k++) {
			if ((U[k + 1] - U[k]) < 0)
			{
				Q[k] = 2 * RO[k] * (U[k + 1] - U[k]) * (U[k + 1] - U[k]);
				//Ð²ÑÐ·ÐºÐ¾ÑÑ‚ÑŒ
			}
			else Q[k] = 0;
		}

		for (int k = 1; k < M - 1; k++)
		{
			U_new[k] = (-tau * (G[k] - G[k - 1]) / deltaS[k]) + U[k];
		}


		//2 Ð·Ð°ÐºÐ¾Ð½
		for (int k = 0; k < M - 1; k++)
		{
			NU_new[k] = NU[k] + tau * (U_new[k + 1] - U_new[k]) / deltaS[k];
			RO_new[k] = 1. / NU[k];
		}
		//Ð²ÑÐ·ÐºÐ¾ÑÑ‚ÑŒ
		for (int k = 0; k < M; k++) {
			if ((U_new[k + 1] - U_new[k]) < 0)
			{
				Q_new[k] = 2 * RO[k] * (U_new[k + 1] - U_new[k]) * (U_new[k + 1] -
					U_new[k]);
			}
			else Q_new[k] = 0;
		}
		//3 Ð·Ð°ÐºÐ¾Ð½
		for (int k = 0; k < M - 1; k++) {
			E_new[k] = ((-tau * (P[k] + 2 * Q[k]) * (((U_new[k + 1] + U[k + 1]) /
				2.0) - ((U_new[k] + U[k]) / 2.0))) / \
				(2 * deltaS[k] + tau * (gamma - 1) * RO_new[k] * (((U_new[k + 1]
					+ U[k + 1]) / 2.0) - ((U_new[k] + U[k]) / 2.0))) + E[k]);
		}
		//Ð£Ð Ð¡
		for (int k = 0; k < M - 1; k++) {
			P_new[k] = (gamma - 1) * E_new[k] / NU_new[k];
		}
		for (int k = 0; k < M - 1; k++) {
			G_new[k] = P_new[k] + Q_new[k];
		}
		double MV = 0, Ek = 0, Ev = 0, E_full = 0, Massa = 0;
		for (int k = 0; k < M; k++)
		{
			MV += U_new[k] * deltaS[k];
			Ek += U_new[k] * U_new[k] * deltaS[k] * 0.5;
			Ev += U_new[k] * E_new[k];
			E_full += Ek + Ev;
			Massa += deltaS[k];
		}
		//fout << Ek <<endl;
		time += tau;
		for (int k = 1; k < M - 1; k++)
		{
			U[k] = U_new[k];
			E[k] = E_new[k];
			NU[k] = NU_new[k];
			G[k] = G_new[k];
			P[k] = P_new[k];

		}
		for (int i = 0; i < M - 1; i++)
		{
			X[i] = X[i] + tau * U_new[i];
		}
		//cout << tau << endl;
	}
	//X_vec.push_back(0);
	//for (int k = 0; k < M - 1; k++)
	//{
	//// deltaX = deltaS / RO_new[k];
	// //X += deltaX;
	// X_vec.push_back(X);
	//}
	//X_vec.push_back(2);
	ofstream fout("test_2.txt");
	//fout << "deltaS: " << deltaS << " tau " << tau << " end_Time " << end_time << endl;
	fout << "X:" << "\t" << "U:" << "\t" << "P:" << "\t" << "RO:" << "\t" << "E:" << "\t" << "Q_new[j]:" << endl;
	int j = 1;
	while (j < M)
	{
		//cout << tau_uv[j]<< endl;
		//fout <<X[j] << endl;
		fout << X[j] << "\t" << U[j] << "\t" << P_new[j] << "\t" << RO_new[j] << "\t" << E[j] << "\t" << Q_new[j] << endl;
		//fout << U[j] << endl;
		j++;
	}
	//fout << X_vec[M] << " " << U[M] <<" " << P[M] <<" " << RO[M] << " "<< E[M] << endl;
	fout.close();
}*/