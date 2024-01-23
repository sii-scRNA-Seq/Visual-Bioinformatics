import { firstValueFrom } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { Injectable } from '@angular/core';

import { Output } from './output';

@Injectable({
  providedIn: 'root'
})
export class UserIdService {

  private userId = '';

  constructor(private http: HttpClient) { 
    this.setUserId();
  }

  async setUserId(): Promise<void> {
    const id = localStorage.getItem('userId') || '';
    const url = 'http://127.0.0.1:5000/getuserid';
    type Params = {[key: string] : string};
    const params: Params = {'user_id' : id};
    const response: Output = await firstValueFrom(this.http.get<Output>(url, {params: params}));
    this.userId = response.text;
    localStorage.setItem('userId', this.userId);
  }
}
