import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { ClientService } from './client.service';
import { UserIdServiceInterface } from './user-id.service.interface';

@Injectable({
  providedIn: 'root'
})
export class UserIdService implements UserIdServiceInterface{
  private readonly userId$: BehaviorSubject<string | null> = new BehaviorSubject<string | null> (null);
  readonly userId: Observable<string | null> = this.userId$.asObservable();

  constructor(private clientService: ClientService) { 
    this.setUserId();
  }

  async setUserId(): Promise<void> {
    const userId = await this.clientService.getUserId(localStorage.getItem('userId') || '');
    localStorage.setItem('userId', userId);
    this.userId$.next(userId);
  }
}