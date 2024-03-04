import { firstValueFrom } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

@Injectable({
  providedIn: 'root'
})
export class BackendHttpClient {

  constructor(private http: HttpClient, private snackBar: MatSnackBar) { }

  async getUserId(id: string): Promise<string> {
    const url = 'http://127.0.0.1:5000/api/getuserid';
    type Params = {[key: string] : string};
    const params: Params = {'user_id': id};
    type Response = {user_id: string};
    let response: Response = {user_id: ''};
    try {
      response = await firstValueFrom(this.http.get<Response>(url, {params: params}));
    } catch (e) {
      this.snackBar.open('Error retrieving userId, please refresh the page and try again', 'Close', { duration: 5000 });
    }
    return response.user_id;
  }
}
