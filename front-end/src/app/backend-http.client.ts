import { firstValueFrom } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { Injectable, isDevMode } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

@Injectable({
  providedIn: 'root'
})
export class BackendHttpClient {

  constructor(private http: HttpClient, private snackBar: MatSnackBar) { }

  async getUserId(id: string): Promise<string> {
    let url = '';
    if (isDevMode()) {
      url = 'http://localhost:5000/api/getuserid';
    } else {
      url = '/api/getuserid';
    }
    type Params = {[key: string] : string};
    const params: Params = {'user_id': id};
    type Response = {user_id: string};
    let response: Response = {user_id: ''};
    try {
      response = await firstValueFrom(this.http.get<Response>(url, {params: params}));
    } catch (e) {
      this.snackBar.open('Error retrieving User ID. Please refresh the page and try again.', 'Close', { duration: 5000 });
    }
    return response.user_id;
  }

  async getDatasetInfo(): Promise<string[]> {
    let url = '';
    if (isDevMode()) {
      url = 'http://localhost:5000/api/getdatasetinfo';
    } else {
      url = '/api/getdatasetinfo';
    }
    type Response = {datasets: string[]};
    let response: Response = {datasets: []};
    try {
      response = await firstValueFrom(this.http.get<Response>(url));
    } catch (e) {
      this.snackBar.open('Error retrieving datasets. Please refresh the page and try again.', 'Close', { duration: 5000 });
    }
    return response.datasets;
  }
}
