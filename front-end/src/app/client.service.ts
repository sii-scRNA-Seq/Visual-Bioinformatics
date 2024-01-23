import { firstValueFrom } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { Injectable } from '@angular/core';

import { Output } from './output';

@Injectable({
  providedIn: 'root'
})
export class ClientService {

  constructor(private http: HttpClient) { }

  async getUserId(id: string): Promise<string> {
    const url = 'http://127.0.0.1:5000/getuserid';
    type Params = {[key: string] : string};
    const params: Params = {'user_id' : id};
    const response: Output = await firstValueFrom(this.http.get<Output>(url, {params: params}));
    return response.text;
  }
}
