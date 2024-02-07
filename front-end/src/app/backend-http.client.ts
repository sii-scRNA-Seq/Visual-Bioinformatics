import { firstValueFrom } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

import { Output } from './output';

@Injectable({
  providedIn: 'root'
})
export class BackendHttpClient {

  constructor(private http: HttpClient, private snackBar: MatSnackBar) { }

  async getUserId(id: string): Promise<string> {
    const url = 'http://127.0.0.1:5000/getuserid';
    type Params = {[key: string] : string};
    const params: Params = {'user_id' : id};
    let response: Output = { text: '' };
    try {
      response = await firstValueFrom(this.http.get<Output>(url, {params: params}));
    } catch (e) {
      this.snackBar.open('Error retrieving userId, please refresh the page and try again', 'Close', { duration: 5000 });
    }
    return response.text;
  }
}
